## Regular Expressions Definition -- regex, regexp
Regular expressions are the patterns used to find or find and replace text on the command line. Regular expressions are used in most modern programming languages and the syntax is usually very similar.

## Regular expressions in Perl

In this tutorial we will be using perl on the command line to showcase how to use regular expressions.  Perl one-liners is a great tool to add to your bioinformatics toolkit to find, replace, extract and generally manipulate strings.

Perl can be executed on the command line like this:
```
perl –ne 'if (/regular_expression/) { print $_.”\n”;}’ input_file_name
```

- `perl` calls the perl program
- `-ne`  these are parameters to the perl program
 - `-n`  feeds the perl program each line of the file
 - `-p`  feeds the perl program each line and prints every line.
 - `-e`  stands for execute what is between the { }
- `if (/regex here/ __)__` is the portion of the perl command where you put a the string you are trying to find.
- `{execute me}`     is the portion of the perl command where you tell perl what to do.  In the example above, each line of the file is printed to the screen.
- `$_` variable in perl that represents the entire line that was fed to the perl program


## Regular Expression Syntax

Like a program language regular expressions have syntax that can represent more than just letters and numbers but specific patterns of letters and numbers.
This syntax will allow you to build a pattern that will match multiple strings or words or sentences.

### Example Text
This file has 12 lines of text containig data about an individuals Name, phone number, address and SS#.
Copy this text into a file named `address.txt`
```
ZEnderX 1-515-999-4321ZZIX
1000 ZBattlesXchool driveZZplayX
SS# 11-27-33-47-57
Bean 1-515-999-3149
1010 Bats drive
SS# 11-22-33-44-55
ZPetraX  1-515-999-1234ZZhouseX
10X0 Battleschool drive ZinX
SS# 21-22-23-24-25
Andrew 1-515-294-1320
206 ZScienceX IZZ.X
SS# 11-02-33-04-50
SS# 11-02-33-04-50
```

### Simple Example that requires no syntax
Recall from above
```
perl -ne 'if (/regular_expression/) { print $_.”\n”;}’ input_file_name
```
To find the line that contains Bean you would type
```
perl -ne 'if (/Bean/) {print $_."\n";}' address.txt
```

Result
```
Bean 1-515-999-3149
```

Other examples
```
perl -ne 'if (/515/) {print $_."\n";}' address.txt
perl -ne 'if (/Battleschool/) {print $_."\n";}' address.txt
perl -ne 'if (/Science/) {print $_."\n";}' address.txt
```

### Example with Simple Syntax `.` , `*`, `?`

the `.` symbol in the regex pattern can represent any character

To find all addresses that have the pattern 10X0 Where X represents any number or character
```
perl -ne 'if(/10.0/) {print $_}' address.txt
```
Result
```
1000 ZBattlesXchool driveZZplayX
1010 Bats drive
10X0 Battleschool drive ZinX
```

To match multiple characters we can use the ```*``` symbol. In the example below it will find all lines with Z followed by 0 or more characters then an X
```
perl -ne 'if(/Z.*X/) {print $_}' address.txt
```
Result prints all lines that contain this pattern 
```
ZEnderX 1-515-999-4321ZIX
1000 ZBattlesXchool driveZZplayX
ZPetraX  1-515-999-1234ZZhouseX
10X0 Battleschool drive ZinX
206 ZScienceX IZZ.X
```

If however we wanted just the characters between Z and X we will need to designate a variable in Perl.  This is achieved by putting the portion of the pattern we want to recall later inside parentheses ```()```
```
perl -ne 'if(/Z(.*)X/) {print $1."\n"}' address.txt
```
Result
```
EnderX 1-515-999-4321ZZI
BattlesXchool driveZZplay
PetraX  1-515-999-1234ZZhouse
in
ScienceX IZZ.
```

Using `.*` will find the pattern with the longest match even if a shorter match is present.  This is called ```Greedy Matching```
Using `.*?` will find the pattern with the shortest match. This is called `Lazy Matching`

```
perl -ne 'if(/Z(.*?)X/) {print $1."\n"}' address.txt
```
Result.  As you can see it finds the shortest match contained between the letters Z and X.
```
Ender
Battles
Petra
in
Science
```

We can add additional variables in the pattern to capture secondary regions.
```
perl -ne 'if(/Z(.*?)X.*ZZ(.*)X/) {print $2."\n"}' address.txt
```
Result.  As you can see this sentence was hard to see until you extracted the text between the pattern ZZ and X.
```
I
play
house
.
```


## Special RegExp Characters

You may want to include these characters in the regular expression you are trying to match.  However, since these characters also have other meanings in perl and programs in general we have to tell the program that this is to be interpreted as a character in the pattern and not part of the syntax of the program.  This is done by "escaping" the special character by preceding it with a back slash (\).  

| Character | Name |
|--------|--------|
| `\` | (backslash) |
| `^` | (caret) |
| `$` | (dollar sign) |
| `.` | (dot) |
| `\|` | (pipe) |
| `?` | (question mark) |
| `*` | (asterisk) |
| `+` | (plus sign) |
| `(` | (open parenthesis) |
| `)` | (closed parenthesis) |
| `[` | (open square bracket) |
| `]` | (closed square bracket) |
| `{` | (open brace) |
| `}` | (closed brace) |


## Character Classes
Below is a table of some of the most common character classes that are used in regular expression.  

| Class | Description |
|--------|------------|
| `[A-Z]` | any single capital letter between A and Z |
| `[a-z]` | any lower case letter between A and Z |
| `[0-9]` | any number between 0 and 9 |
| `\s` | matches any whitespace character |
| `\w` | matches a word character equivelent to [A-Za-z0-9_] |
| `[:ascii:]` | Matches any character in the ASCII character set |
| `\d` | matches any number equivalent to [0-9] |
| `\d{3}` | matches any string of 3 numbers |
| `\d{2,5}` | matches any string of 2 to 5 numbers |

A more comprehensive list of character classes can be found here.

[http://perldoc.perl.org/perlrecharclass.html](http://perldoc.perl.org/perlrecharclass.html)

