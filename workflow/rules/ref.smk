rule bwameth_index:
    input:
        ancient(REFERENCE_GENOME)
    output:
        REFERENCE_GENOME + '.bwameth.c2t',
        REFERENCE_GENOME + '.bwameth.c2t.amb',
        REFERENCE_GENOME + '.bwameth.c2t.ann',
        REFERENCE_GENOME + '.bwameth.c2t.bwt',
        REFERENCE_GENOME + '.bwameth.c2t.pac',
        REFERENCE_GENOME + '.bwameth.c2t.sa'
    params:
        bwameth_path = config['paths']['bwameth_path'],
    shell:
        'bwameth.py index {input}'