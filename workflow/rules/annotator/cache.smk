rule download_snpeff_cache:
    conda:
        "../../envs/awscli.yaml"
    output:
        dir=directory(path_cache_snpeff),
    params:
        dir_src=f"{config['genome']}.{config['version_snpeff']}",
        dir_trgt=config["cache_snpeff"],
    log:
        "logs/download_snpeff_cache.log",
    shell:
        """
        {{ url_base="s3://annotation-cache/snpeff_cache"
        caches_aws=( $(aws s3 ls --no-sign-request ${{url_base}}/ | grep PRE | sed 's/^[[:space:]]*//' | cut -d ' ' -f 2 | sed 's/\\/$//' | tr '\n' ' ') )

        if [[ " ${{caches_aws[@]}} " =~ " {params.dir_src} " ]]; then
            aws s3 --no-sign-request sync ${{url_base}}/{params.dir_src} {params.dir_trgt}
        else
            snpEff download {params.dir_src} -dataDir {params.dir_trgt}
        fi; }} \\
        1> {log} 2>&1
        """


rule download_vep_cache:
    conda:
        "../../envs/awscli.yaml"
    output:
        dir=directory(path_cache_vep),
    params:
        dir_src=f"{config['version_vep']}_{config['genome']}",
        dir_trgt=f"{config['cache_vep']}/{config['species']}",
        cache=config["cache_vep"],
        species=config["species"],
        version=config["version_vep"],
        genome=config["genome"],
    log:
        "logs/download_vep_cache.log",
    shell:
        """
        {{ url_base="s3://annotation-cache/vep_cache"
        caches_aws=( $(aws s3 ls --no-sign-request ${{url_base}}/ | grep PRE | sed 's/^[[:space:]]*//' | cut -d ' ' -f 2 | sed 's/\\/$//' | tr '\n' ' ') )

        if [[ " ${{caches_aws[@]}} " =~ " {params.dir_src} " ]]; then
            aws s3 --no-sign-request sync ${{url_base}}/{params.dir_src} {params.dir_trgt}
        else
            vep_install \\
            --CACHEDIR {params.cache} \\
            --DESTDIR {params.cache} \\
            --PLUGINSDIR {params.cache}/Plugins/ \\
            --CACHE_VERSION {params.version} \\
            --SPECIES {params.species} \\
            --ASSEMBLY {params.genome} \\
            --PREFER_BIN --NO_UPDATE --AUTO cf
        fi; }} \\
        > {log} 2>&1
        """


rule download_annotsv_cache:
    conda:
        "../../envs/git.yaml"
    output:
        **get_annotsv_cache_outputs(),
        annotsv=temp(directory(f"{config['cache_annotsv']}/AnnotSV")),
    params:
        **get_annotsv_cache_parameters(),
        dir=config["cache_annotsv"],
        version=config["version_annotsv"],
    log:
        "logs/download_annotsv_cache.log",
    shell:
        """
        {{ git clone https://github.com/lgmgeo/AnnotSV.git {output.annotsv}
        cd {output.annotsv} || exit 1
        git checkout {params.version}

        make PREFIX=. install
        make PREFIX=. {params.arg_install}

        for dir in {params.dirs}; do
            mv ${{dir}} ../
        done; }} \\
        > {log} 2>&1
        """


rule query_annovar_protocols:
    output:
        txt=f"{config['cache_annovar']}/{GENOME2}_avdblist.txt",
    log:
        "logs/query_annovar_protocols.log",
    params:
        dir=config["cache_annovar"],
        genome=GENOME2,
    shell:
        """
        {{ annotate_variation.pl \\
            -buildver {params.genome} \\
            -webfrom annovar \\
            -downdb avdblist \\
            {params.dir} }} \\
        > {output.txt} 2> {log}
        """


rule download_annovar_cache:
    input:
        txt=ancient(f"{config['cache_annovar']}/{GENOME2}_avdblist.txt"),
    output:
        txt=f"{config['cache_annovar']}/{GENOME2}_{{protocol}}.txt",
    params:
        dir=config["cache_annovar"],
        file=f"{GENOME2}_{{protocol}}.txt.gz",
        genome=GENOME2,
        is_annovar=lambda wildcards: wildcards.protocol not in PROTOCOLS_UCSC,
        arg_webfrom=lambda wildcards: (
            f"-webfrom annovar" if wildcards.protocol not in PROTOCOLS_UCSC else ""
        ),
    log:
        "logs/download_annovar_cache.{protocol}.log",
    run:
        if params.is_annovar:
            with open(input.txt) as f:
                avdblist = {line.split()[0] for line in f if line.strip()}

            if params.file not in avdblist:
                raise ValueError(
                    f"Protocol {wildcards.protocol} not found in avdblist."
                )

        shell(
            """
            annotate_variation.pl \\
                -buildver {params.genome} \\
                -downdb \\
                {params.arg_webfrom} \\
                {wildcards.protocol} \\
                {params.dir} \\
            > {log} 2>&1
            """
        )
