
FN=marathon

pdf:
	pandoc -o $(FN).pdf $(FN).md --latex-engine=xelatex --template=templates/marathon.latex -V geometry:margin=1.2in

docx:
	pandoc -o $(FN).docx $(FN).md
