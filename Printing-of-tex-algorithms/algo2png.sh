#!/bin/bash
# algo2png.sh
BASENAME=$(basename "$1" .tex)

cat > _temp.tex << 'EOF'
\documentclass{article}
\usepackage{algorithm,algorithmic,amsmath,amssymb,mathrsfs}
\usepackage[active,tightpage]{preview}
\PreviewEnvironment{algorithm}
\begin{document}
EOF
cat "$1" >> _temp.tex
cat >> _temp.tex << 'EOF'
\end{document}
EOF
pdflatex -interaction=nonstopmode _temp.tex > /dev/null
pdftoppm -png -singlefile _temp.pdf "$BASENAME"
rm _temp.*