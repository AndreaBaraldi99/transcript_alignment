### Program that counts how many reads support a transcript

This program takes a bam file with some reads and a gtf file with transcripts and exons and checks whether the reads support the transcripts.
It does so by checking if the exons and introns regions of the reads align correctly to the regions of the transcript (with a 1 base editable tollerance).
At the end, it prints out the count of the reads that support a transcript for every transcript inside of the gtf file
