# This Makefile lazily/simply uses latexmk!

THESIS = thesis

.PHONY: FORCE_MAKE draft final
${THESIS}.pdf: ${THESIS}.tex FORCE_MAKE
	latexmk --pdf $<

draft: ${THESIS}.pdf

final:
	latexmk --pdf --jobname=${THESIS}-final \
		'\\def\\FinalResults{1}\\input{${THESIS}.tex}'

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
