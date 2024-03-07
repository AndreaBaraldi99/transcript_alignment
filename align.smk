rule count_reads:
    input:
        bam_file='sim.ENST00000466434.Aligned.sortedByCoord.out.bam',
        gtf_file='ENST00000466434.gtf'
    params:
        quality = 0
    shell:
        '''
        python3 align.py --bam {input.bam_file} --gtf {input.gtf_file} --quality {params.quality}
        '''
