%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% homus_workus
% LaTeX Template
% Version 1.1 (September 14, 2024)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ------------------------------------------------------------
%	Class initialisation
% ------------------------------------------------------------
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{HomusWorkus}[2024/09/14 Example LaTeX class]
\LoadClass{article}
% ------------------------------------------------------------
%                       Size of paper
% ------------------------------------------------------------
\usepackage[
    a4paper,
    top=2cm,
    bottom=2cm,
    left=3cm,
    right=3cm,
    marginparwidth=1.75cm
]{geometry}
% ------------------------------------------------------------
%                       Language and font settings
% ------------------------------------------------------------
\usepackage[T2A,T2A]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[english,ukrainian]{babel}

% \babelfont{rm}{FreeSerif}
% \babelfont{sf}{FreeSans}
% \babelfont{tt}{FreeMono}
% ------------------------------------------------------------
%                       Packages for contents
% ------------------------------------------------------------
\usepackage[colorlinks=true, allcolors=black]{hyperref} % table of contsnts
% ------------------------------------------------------------
%                       Math settings
% ------------------------------------------------------------
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage{algorithm2e}
\usepackage{listings}
% ------------------------------------------------------------
%                       Define colours
% ------------------------------------------------------------
\RequirePackage{xcolor}
\definecolor{slcolor}{HTML}{882B21}
\definecolor{main}{HTML}{5989cf}    % setting main color to be used
\definecolor{sub}{HTML}{cde4ff}     % setting sub color to be used
% ------------------------------------------------------------
%                       Graphics settings
% ------------------------------------------------------------
\usepackage{graphicx}
\usepackage{tikz}
\usetikzlibrary{trees}
%\usetikzlibrary{shapes,arrows,positioning}
\usetikzlibrary{automata,arrows,positioning,calc}
% ------------------------------------------------------------
\usepackage{multicol}       % Fro multi column page
\usepackage{tabularx}
\usepackage{enumitem}       % set lists properties 
\setlist[itemize]{noitemsep, topsep=0pt}
\setlist[enumerate]{noitemsep, topsep=0pt}
% ------------------------------------------------------------
%                       Links
% ------------------------------------------------------------
\labelformat{algocf}{(\textit{alg.}\,#1)}
\labelformat{table}{(\textit{tab.}\,#1)}
\labelformat{figure}{(\textit{fig.}\,#1)}
% ------------------------------------------------------------
%                       End
% ------------------------------------------------------------

