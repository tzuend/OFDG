THESIS = thesis
LATEX = pdflatex -interaction=nonstopmode -halt-on-error
DRAFT_PRODUCTS = ${THESIS}.aux ${THESIS}.bbl ${THESIS}.blg \
	${THESIS}.fdb_latexmk ${THESIS}.fls ${THESIS}.log ${THESIS}.out \
	${THESIS}.pdf ${THESIS}.thm ${THESIS}.toc
FINAL_PRODUCTS = ${THESIS}-final.aux ${THESIS}-final.bbl ${THESIS}-final.blg \
	${THESIS}-final.fdb_latexmk ${THESIS}-final.fls ${THESIS}-final.log \
	${THESIS}-final.out ${THESIS}-final.pdf ${THESIS}-final.thm ${THESIS}-final.toc

.PHONY: draft final verify-draft verify-final

draft:
	rm -f ${DRAFT_PRODUCTS}
	${LATEX} ${THESIS}.tex
	bibtex ${THESIS}
	${LATEX} ${THESIS}.tex
	${LATEX} ${THESIS}.tex
	${LATEX} ${THESIS}.tex
	$(MAKE) verify-draft

final:
	rm -f ${FINAL_PRODUCTS}
	${LATEX} --jobname=${THESIS}-final \
		'\def\FinalResults{1}\input{${THESIS}.tex}'
	bibtex ${THESIS}-final
	${LATEX} --jobname=${THESIS}-final \
		'\def\FinalResults{1}\input{${THESIS}.tex}'
	${LATEX} --jobname=${THESIS}-final \
		'\def\FinalResults{1}\input{${THESIS}.tex}'
	$(MAKE) verify-final

verify-draft:
	@if grep -Eq 'LaTeX Warning: (Citation|Reference).*undefined|LaTeX Warning: There were undefined references|Overfull \\[hv]box' ${THESIS}.log; then \
		echo "Draft report contains unresolved references or overfull boxes."; \
		exit 1; \
	fi

verify-final:
	@if grep -Eq 'LaTeX Warning: (Citation|Reference).*undefined|LaTeX Warning: There were undefined references|Overfull \\[hv]box' ${THESIS}-final.log; then \
		echo "Final report contains unresolved references or overfull boxes."; \
		exit 1; \
	fi

.PHONY: clean
clean:
	latexmk -c -bibtex ${THESIS}

cleanall: clean
	latexmk -C -bibtex ${THESIS}

# Only useful to create a template!
thesis-template.zip: Makefile latexmkrc *.tex ETHlogo.* refs.bib \
                     aas_macros.sty declaration-originality.pdf
	mkdir -p thesis
	cp -v Makefile latexmkrc *.tex ETHlogo.* refs.bib aas_macros.sty \
              declaration-originality.pdf thesis
	cp -v declaration-originality.pdf                 \
              thesis/declaration-originality_original.pdf
	rm -f thesis-template.zip
	zip -r thesis-template.zip thesis/
	rm -r thesis
