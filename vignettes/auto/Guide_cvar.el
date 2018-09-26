(TeX-add-style-hook
 "Guide_cvar"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("jss" "nojss" "article")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("geometry" "left=2cm" "right=2cm" "bottom=15mm") ("natbib" "authoryear" "round" "longnamesfirst")))
   (add-to-list 'LaTeX-verbatim-environments-local "alltt")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "jss"
    "jss10"
    "fontenc"
    "geometry"
    "graphicx"
    "color"
    "alltt"
    "natbib"
    "hyperref"
    "amsmath"
    "amsfonts")
   (TeX-add-symbols
    '("ES" ["argument"] 1)
    '("VaR" ["argument"] 1))
   (LaTeX-add-labels
    "sec:introduction"
    "sec:var"
    "sec:expected-shortfall")
   (LaTeX-add-bibliographies
    "REFERENCES"))
 :latex)

