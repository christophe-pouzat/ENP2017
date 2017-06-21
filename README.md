# Course material for the ENP June 20 2017 lecture: "Fluorescence Imaging Analysis: The Case of Calcium Transients."

The slides of the lecture are in file `Pouzat_ENP2017.pdf`. The handout, nicer for printing, with some additional stuff is in file `Pouzat_ENP2017_handout.pdf`. All the analysis is done with [`C`](https://en.wikipedia.org/wiki/C_(programming_language)) together with the [GSL](https://www.gnu.org/software/gsl/) library (the Gnu Scientific Library), while figures are done with [gnuplot](http://gnuplot.info/). The `C` code together with the full details of the whole analysis can be found in file `ENP2017_technical_notes.html`. The `C` source files are in directory `code`. All figures except one are created as part of the analysis--the `gnuplot` code generating the figures is therefore found in `ENP2017_technical_notes.html`--.

Each of the three file mentioned above are generated from their [org](http://orgmode.org/) source file (files with the same prefix terminating with the `.org` suffix). The `C` code is presented following Don Knuth's [literate programming](https://en.wikipedia.org/wiki/Literate_programming) paradigm implemented in `org` (see the paper of Eric Schulte, Dan Davison, Thomas Dye, Carsten Dominik [A Multi-Language Computing Environment for Literate Programming and Reproducible Research](https://www.jstatsoft.org/article/view/v046i03), 2012, to understand what that means and how that works).

The only actual code contained in `Pouzat_ENP2017.org` and `Pouzat_ENP2017_handout.pdf` shows how to automatically download a `pdf` file from the web, open it at a specific page and extract a specific figure (more precisely a specific part of the page)--all that exclusively from the command line--. This is used to get the only figure of the lecture that is not a result of the analysis.
