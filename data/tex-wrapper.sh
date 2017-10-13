#! /bin/bash
# This script simply adds wrapping lines to the plot tex files to enable standalone compilation and then compiles it

# Prepend lines using sed. Source: https://www.cyberciti.biz/faq/bash-prepend-text-lines-to-file/
sed -i '1s;^;\\\documentclass{standalone}\n\\usepackage{tikz}\n\\\begin{document}\n;' $1

# Append line at end to close document
sed -i '$ a \\\end{document}' $1

# Compile the file
pdflatex $1

# Construct pdf file name from tex file name and show it
pdf=$(echo $1 | sed 's/\..*$/\.pdf/')
vupdf $pdf
