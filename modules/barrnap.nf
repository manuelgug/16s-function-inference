// rRNA prediction (barrnap 0.9), one genome per CPU. Kingdom from the GTDB domain (column 3).
process BARRNAP {
    tag "${batch}"

    input:
    tuple val(batch), path(genomes), path(fna)

    output:
    tuple val(batch), path('gff')

    script:
    """
    mkdir gff
    tail -n +2 ${genomes} | cut -f1,3 | while read -r acc domain; do
        [ -s "${fna}/\$acc.fna.gz" ] || continue
        if [ "\$domain" = Archaea ]; then echo "\$acc arc"; else echo "\$acc bac"; fi
    done | xargs -r -P ${task.cpus} -n 2 sh -c \\
        'gunzip -c ${fna}/\$0.fna.gz > \$0.fna && barrnap --quiet --kingdom \$1 --threads 1 \$0.fna > gff/\$0.gff && rm -f \$0.fna \$0.fna.fai'
    """
}
