
## Learning Objective
Upon completion of this section on fastq quality scores the learner will understand the following:
* Assci character encoding are used to represent numbers
* These numbers are converted to values between 0 and 41 to represent quality score


<pre>
  <span style="color: purple">SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS</span>.....................................................
  ..........................<span style="color: green">XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX</span>......................
  ...............................<span style="color: blue">IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII</span>......................
  .................................<span style="color: orange"><b>J</b>JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ</span>......................
  <span style="color: brown">LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL</span>....................................................
  !"#$%&amp;'()*+,-./0123456789:;&lt;=&gt;?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
  |                         |    |        |                              |                     |
 33                        59   64       73                            104                   126
<span style="color: purple">  0........................26...31.......40                                </span>
<span style="color: green">                           -5....0........9.............................40 </span>
<span style="color: blue">                                 0........9.............................40 </span>
<span style="color: orange">                                    3.....9.............................40 </span>
<span style="color: brown">  0.2......................26...31........41                              </span>

 <span style="color: purple">S - Sanger        Phred+33,  raw reads typically (0, 40)</span>
 <span style="color: green">X - Solexa        Solexa+64, raw reads typically (-5, 40)</span>
 <span style="color: blue">I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)</span>
 <span style="color: orange">J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
     with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
     (Note: See discussion above).</span>
 <span style="color: brown">L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)</span>
</pre>


[This table was taken from ](https://en.wikipedia.org/wiki/FASTQ_format)where more information can be found on this topic.

To determine if they score is <blockcode>Phred+33</blockcode>, <blockcode>Phred+64</blockcode> or <blockcode>Solexa+64</blockcode>, use this one-liner (you can use <blockcode>zcat</blockcode>  if the <blockcode>fastq</blockcode> file is gzipped):

```
head -n 10000 input.fastq |\
  awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | \
  awk 'BEGIN{min=100;max=0;} \
      {for(i=1;i<=NF;i++) \
          {if($i>max) max=$i; \
               if($i<min) min=$i;}}END \
          {if(max<=74 && min<59) \
                     print "Phred+33"; \
           else \
           if(max>73 && min>=64) \
                     print "Phred+64"; \
           else \
           if(min>=59 && min<64 && max>73) \
                     print "Solexa+64"; else print "Unknown score encoding!";}'
 ```

 [Table of contents](/index.md)
