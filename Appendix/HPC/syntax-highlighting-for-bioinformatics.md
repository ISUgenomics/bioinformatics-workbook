# Syntax high-lighting for Biology specific files

If you're a coder, then you already know how useful the syntax highlighting is for your programming language. However, if you're a biologist and work with lots of biology specific files (fasta, msa, vcf etc) on command-line, then you might have noticed how non-intuitive it feels for manually inspecting them. This tutorial aims to reduce that feeling a little bit and makes working on CLI a bit better! If you're a fan of `nano` for editing files, check out the first section. If you use `vim` then move to the second section. For now, we are leaving out the `emacs` but will be adding eventually.

## Nano as the text editor

For `nano` most settings can be set in `.nanorc` file located in your home directory (`/home/username/.nanorc` or simply `~/.nanorc`).

You can add file specific syntax coloring by editing `.nanorc` file. For example, for nucleotide file (file extensions `.fna` `.fasta` or `.fa`), you can add this section:

```bash
#colored based on purines or pyrimidines
syntax "fasta" "\.fna$" "\.fasta$" "\.fa$"
color brightblue "[AT]"
color brightgreen "[GC]"
color brightred "[N]"
```


Similarly, for amino acid fasta files, you can add these lines (ending with `.faa` or `.aa`):


```bash
#colored based on aminoacid properties
syntax "fasta" "\.faa$" "\.aa$"
color brightblue "[AILMFWV]"
color brightred "[RK]"
color brightgreen "[NQ]"
color white "[C]"
color magenta "[ED]"
color red "[G]"
color cyan "[HY]"
color brightyellow "[P]"
color green "[ST]"
```

This can also be used for alignment files (nucleotide or proteins) by adding more extensions to the above syntax coloring.
