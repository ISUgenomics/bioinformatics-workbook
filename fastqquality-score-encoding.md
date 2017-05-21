
{{ style("color: purple") "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS}}.....................................................
..........................{{ style("color: green") "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX}}......................
...............................{{ style("color: blue") "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII}}......................
.................................{{ style("color: orange") "<b>J</b>JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ}}......................
{{ style("color: brown") "LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL}}....................................................
!"#$%&amp;'()*+,-./0123456789:;&lt;=&gt;?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
| | | | | |
33 59 64 73 104 126
{{ "0........................26...31.......40" | style(“color: purple") }} }}
{{ "-5....0........9.............................40” | style("color: green")  }}
{{ "0........9.............................40” | style("color: blue")  }}
{{ "3.....9.............................40” | style("color: orange")  }}
{{ "0.2......................26...31........41” | style("color: brown")  }}

{{ "S - Sanger Phred+33, raw reads typically (0, 40)” | style("color: purple") }}
{{ “X - Solexa Solexa+64, raw reads typically (-5, 40)” | style("color: green") }}
{{ "I - Illumina 1.3+ Phred+64, raw reads typically (0, 40)” | style("color: blue") }}
{{ "J - Illumina 1.5+ Phred+64, raw reads typically (3, 40)
with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
(Note: See discussion above).” | style("color: orange") }}
{{ "L - Illumina 1.8+ Phred+33, raw reads typically (0, 41)” | style("color: brown") }}


[For more information see Wikipedia entry](https://en.wikipedia.org/wiki/FASTQ_format)
