(TeX-add-style-hook "CrossScreening-vignette"
 (lambda ()
    (LaTeX-add-bibliographies
     "ref")
    (LaTeX-add-labels
     "sec:intro"
     "sec:obs-study"
     "sec:sen"
     "sec:sen-value"
     "sec:appr-power-sens"
     "sec:cross-screen"
     "sec:discussion")
    (TeX-run-style-hooks
     "cleveref"
     "authblk"
     "natbib"
     "authoryear"
     "amsmath"
     "url"
     "geometry"
     ""
     "fontenc"
     "T1"
     "mathpazo"
     "sc"
     "latex2e"
     "art11"
     "article"
     "11pt")))

