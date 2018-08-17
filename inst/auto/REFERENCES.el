(TeX-add-style-hook
 "REFERENCES"
 (lambda ()
   (LaTeX-add-bibitems
    "acerbi2002expected"
    "PerformanceAnalytics2018"
    "VaRES2013"
    "actuarJSS2008")
   (LaTeX-add-environments
    '("eptblFigure" LaTeX-env-args ["argument"] 0)
    '("epfigFigure" LaTeX-env-args ["argument"] 0)
    '("epFigure" LaTeX-env-args ["argument"] 0)))
 :bibtex)

