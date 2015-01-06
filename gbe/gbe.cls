%%
%% This is file `gbe.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% classes.dtx  (with options: `article')
%%
%% This is a generated file.
%%
%% Copyright 1993 1994 1995 1996 1997 1998 1999 2000 2001
%% The LaTeX3 Project and any individual authors listed elsewhere
%% in this file.
%%
%% This file was generated from file(s) of the LaTeX base system.
%% --------------------------------------------------------------
%%
%% It may be distributed and/or modified under the
%% conditions of the LaTeX Project Public License, either version 1.2
%% of this license or (at your option) any later version.
%% The latest version of this license is in
%%    http://www.latex-project.org/lppl.txt
%% and version 1.2 or later is part of all distributions of LaTeX
%% version 1999/12/01 or later.
%%
%% This file may only be distributed together with a copy of the LaTeX
%% base system. You may however distribute the LaTeX base system without
%% such generated files.
%%
%% The list of all files belonging to the LaTeX base distribution is
%% given in the file `manifest.txt'. See also `legal.txt' for additional
%% information.
%%
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{gbe}
              [2008/07/24 v1.4e
 Standard LaTeX document class]
\newcommand\@ptsize{}
\newif\if@restonecol
\newif\if@titlepage
\@titlepagefalse
\if@compatibility\else
\DeclareOption{a4paper}
   {\setlength\paperheight {297mm}%
    \setlength\paperwidth  {210mm}}
\DeclareOption{a5paper}
   {\setlength\paperheight {210mm}%
    \setlength\paperwidth  {148mm}}
\DeclareOption{b5paper}
   {\setlength\paperheight {250mm}%
    \setlength\paperwidth  {176mm}}
\DeclareOption{letterpaper}
   {\setlength\paperheight {11truein}%
    \setlength\paperwidth  {8.5truein}}
\DeclareOption{legalpaper}
   {\setlength\paperheight {14in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{executivepaper}
   {\setlength\paperheight {10.5in}%
    \setlength\paperwidth  {7.25in}}
\DeclareOption{landscape}
   {\setlength\@tempdima   {\paperheight}%
    \setlength\paperheight {\paperwidth}%
    \setlength\paperwidth  {\@tempdima}}
\fi
\if@compatibility
  \renewcommand\@ptsize{0}
\else
\DeclareOption{10pt}{\renewcommand\@ptsize{0}}
\fi
\DeclareOption{11pt}{\renewcommand\@ptsize{1}}
\DeclareOption{12pt}{\renewcommand\@ptsize{2}}
\if@compatibility\else
\DeclareOption{oneside}{\@twosidefalse \@mparswitchfalse}
\fi
\DeclareOption{twoside}{\@twosidetrue  \@mparswitchtrue}
\DeclareOption{draft}{\setlength\overfullrule{1pt}}
\if@compatibility\else
\DeclareOption{final}{\setlength\overfullrule{0pt}}
\fi
\DeclareOption{titlepage}{\@titlepagetrue}
\if@compatibility\else
\DeclareOption{notitlepage}{\@titlepagefalse}
\fi
\if@compatibility\else
\DeclareOption{onecolumn}{\@twocolumnfalse}
\fi
\DeclareOption{twocolumn}{\@twocolumntrue}
\DeclareOption{leqno}{\input{leqno.clo}}
\DeclareOption{fleqn}{\input{fleqn.clo}}
\DeclareOption{openbib}{%
  \AtEndOfPackage{%
   \renewcommand\@openbib@code{%
     \advance\leftmargin\bibindent
      \itemindent -\bibindent
      \listparindent \itemindent
      \parsep \z@
      }%
   \renewcommand\newblock{\par}}%
}



\def\igDelta#1{
\let\symb=#1
}

%\DeclareOption{mathptmx}{\AtEndOfClass{\RequirePackage[cmbold,subscriptcorrection,nofontinfo]{mathptmx}}

%\DeclareMathSymbol\Gamma  {\mathord}{letters}{48}
%\ifx Y\symb
%\DeclareMathSymbol\Delta  {\mathord}{letters}{49}
%\fi
%\DeclareMathSymbol\Theta  {\mathord}{letters}{50}
%\DeclareMathSymbol\Lambda {\mathord}{letters}{51}
%\DeclareMathSymbol\Xi     {\mathord}{letters}{52}
%\DeclareMathSymbol\Pi     {\mathord}{letters}{53}
%\DeclareMathSymbol\Sigma  {\mathord}{letters}{54}
%\DeclareMathSymbol\Upsilon{\mathord}{letters}{55}
%\DeclareMathSymbol\Phi    {\mathord}{letters}{56}
%\DeclareMathSymbol\Psi    {\mathord}{letters}{57}
%\DeclareMathSymbol\Omega  {\mathord}{letters}{127}
%}
\DeclareOption{cmfonts}{%
    \g@addto@macro{\@greeksetup}{%
        \DeclareMathSymbol\varGamma    {\mathord}{letters}{0}
    \ifx Y\symb
       \DeclareMathSymbol\varDelta    {\mathord}{letters}{1}
     \fi
          \DeclareMathSymbol\varTheta    {\mathord}{letters}{2}
       \DeclareMathSymbol\varLambda   {\mathord}{letters}{3}
       \DeclareMathSymbol\varXi       {\mathord}{letters}{4}
       \DeclareMathSymbol\varPi       {\mathord}{letters}{5}
       \DeclareMathSymbol\varSigma    {\mathord}{letters}{6}
       \DeclareMathSymbol\varUpsilon  {\mathord}{letters}{7}
       \DeclareMathSymbol\varPhi      {\mathord}{letters}{8}
       \DeclareMathSymbol\varPsi      {\mathord}{letters}{9}
       \DeclareMathSymbol\varOmega    {\mathord}{letters}{10}
      }
}
\DeclareOption{slantedgreek}{\def\@greeksetup{\AtBeginDocument{%
\ifx Y\symb
\renewcommand{\Gamma}{{\mit\varGamma}}

\renewcommand{\Delta}{{\mit\varDelta}}

\renewcommand{\Theta}{{\mit\varTheta}}
\renewcommand{\Lambda}{{\mit\varLambda}}
\renewcommand{\Xi}{{\mit\varXi}}
\renewcommand{\Pi}{{\mit\varPi}}
\renewcommand{\Sigma}{{\mit\varSigma}}
\renewcommand{\Upsilon}{{\mit\varUpsilon}}
\renewcommand{\Phi}{{\mit\varPhi}}
\renewcommand{\Psi}{{\mit\varPsi}}
\renewcommand{\Omega}{{\mit\varOmega}}
\fi}}}
\igDelta{Y}

\DeclareOption{uprightgreek}{\def\@greeksetup{}}

\ExecuteOptions{letterpaper,10pt,twoside,twocolumn,final,uprightgreek,cmfonts}

\ProcessOptions
\input{size1\@ptsize.clo}

\usepackage{lastpage,fancybox,}


%%%%%%%%%%%%%%%%%%%%trim size setting

%xxx sdj \textwidth497.9232pt%
%\textheight626.1454pt
%xxx sdj \textheight646.1454pt
\columnwidth238pt
\columnsep22.7622pt
%xxx sdj \oddsidemargin81pt%
%xxx sdj \evensidemargin78pt%5mm
%xxx sdj \topmargin42pt%
%xxx sdj \headheight10pt
\headsep 36pt %xxx sdj \headsep75.4pt%
\parindent11pt
\footskip 44pt%
\setlength{\fboxsep}{0pt}
\@greeksetup

%%%%%%%%%%%%%%%Plain Style


\def\label@firstpage{%
   \immediate\write\@auxout{\string
   \newlabel{firstpage}{{}{\thepage}}}}
  \newlabel{firstpage}{{}{1}}

\def\public#1{\let\@pg=#1}
\def\pageproof#1{\let\@pp=#1}
\def\jno#1{\gdef\@kxi{#1}}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\ps@plaingbe{%
     \def\@oddhead{\CROPMARKSA\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lr@{}}\\[3pt]
     \hfill  \lrjfont \textcolor{runheadcolor}{GBE}\\[-1.5pt]
     \hskip-0pt\textcolor{sectioncolor}{\rule{\textwidth}{.3pt}} \end{tabular*}}%


%xxx sdj \def\@oddfoot{\gridpl\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lr@{}}\\[-40pt]
%xxx sdj {\copyfont{$\gbecopyright$} The Author(s) 2010. Published by Oxford University
\def\@oddfoot{\gridpl\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lr@{}}\\
{\copyfont{\copyright} The Author(s) 2010. Published by Oxford University
Press on behalf of the \textit{Society for Molecular Biology and Evolution.}}\\[-2.7pt]
{\copyfont This\hspace*{0.5pt} is\hspace*{0.5pt} an\hspace*{0.5pt} Open\hspace*{0.3pt} Access\hspace*{0.5pt} article\hspace*{0.5pt} distributed\hspace*{0.5pt} under\hspace*{0.3pt} the\hspace*{0.5pt} terms\hspace*{0.5pt} of\hspace*{0.5pt}
the\hspace*{0.5pt} Creative\hspace*{0.5pt} Commons\hspace*{0.5pt} Attribution\hspace*{0.3pt} Non-Commercial\hspace*{0.5pt} License\hspace*{0.5pt} (http://creativecommons.org/licenses/by-nc/}\\[-2.7pt]
{\copyfont 2.5), which permits unrestricted non-commercial use, distribution,
and reproduction in any medium, provided the original work is properly cited.}\hfill\\[8.5pt]
\textcolor{sectioncolor}{\rule{\textwidth}{0.3pt}}\\[4pt]
\catchfont{\textit{Genome Biol. Evol.}}\hspace*{6.4pt}\@vol:\fpage--\lpage.\quad\@rundoi\quad\@advance\hfill
\firstfont{\thepage}\vspace*{-4pt}\\[86pt]\end{tabular*}}

        \def\@evnhead{\CROPMARKSA}%
      \ifx N\@pp%
     \def\@oddhead{\CROPMARKSA\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lr@{}}\\[3pt]
\hfill  \lrjfont \textcolor{runheadcolor}{GBE}\\[-1.5pt]
\hskip-0pt\textcolor{sectioncolor}{\rule{\textwidth}{.3pt}} \end{tabular*}}%

      \else
      \ifx N\@pg


     \def\@oddhead{\CROPMARKSA\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lr@{}}\\[3pt]
     \hfill  \lrjfont \textcolor{runheadcolor}{GBE}\\[-1.5pt]
     \hskip-0pt\textcolor{sectioncolor}{\rule{\textwidth}{.3pt}} \end{tabular*}}%

      \else

     \def\@oddhead{\CROPMARKSA\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lr@{}}\\[3pt]
\hfill  \lrjfont \textcolor{runheadcolor}{GBE}\\[-1.5pt]
\hskip-0pt\textcolor{sectioncolor}{\rule{\textwidth}{.3pt}} \end{tabular*}}%
      \fi
      \fi
       \let\@mkboth\markboth
}




\newcommand{\journal}[6]{\def\@jnone{#1}\def\@jntwo{#2}\def\@jnthree{#3}\def\@jnfour{#4}\def\@jnfive{#5}\def\@jnsix{#6}}
\newcommand{\copyeditor}[1]{\def\@ce{#1}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\def\keywords#1{{\def\thefootnote{}\footnotetext{\vspace*{3.5pt}\hskip 1pt {Key words}: #1}}}
\def\email#1{{\def\thefootnote{}\footnotetext{\hskip 1pt E-mail: #1}}}
\def\jourinfo{{\def\thefootnote{}\footnotetext{\\[-5.57pt]
\hskip -14pt\fontsize{8.3}{11.2}\norfont  {\it Mol. Biol. Evol}. \@vol(\@iss):\fpage{}--\lpage. 2009\\
        doi:10.1093/molbev/\@kxi\\
        \@advance} }}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%grid for textwidth 140mm
\def\gridframe#1{
\let\grid=#1
\ifx Y\grid
\def\gridpl{\setlength{\fboxsep}{0pt}\setlength{\unitlength}{1pc}\begin{picture}(0,0)
\put(0,5.35){\framebox(41.5,53.55)}%33.1,47
\end{picture}}

%\def\gridpl{\setlength{\fboxsep}{0pt}\setlength{\unitlength}{1pc}\begin{picture}(0,0)
%\put(0,2.3){\framebox(41.5,63.1)}%33.1,47
%\end{picture}}


\def\gridod{\setlength{\unitlength}{1pc}\begin{picture}(0,0)
\put(0,-58.9){\framebox(41.5,53.55)}%-46.458
\end{picture}}

%\def\gridod{\setlength{\unitlength}{1pc}\begin{picture}(0,0)
%\put(0,-61.95){\framebox(41.5,63.2)}%-46.458
%\end{picture}}


\def\gridev{\setlength{\unitlength}{1pc}\begin{picture}(0,0)
\put(0,-58.9){\framebox(41.5,53.55)}%-46.458
\end{picture}


%\def\gridev{\setlength{\unitlength}{1pc}\begin{picture}(0,0)
%\put(0,-61.95){\framebox(41.5,63.2)}%-46.458
%\end{picture}

}
\else
\def\gridpl{}
\def\gridod{}
\def\gridev{}
\fi
}


\renewcommand\normalsize{%
  \@setfontsize\normalsize{9.5}{12}\norfont%
  \abovedisplayskip 6\p@ \@plus 2\p@ \@minus5\p@
  \abovedisplayshortskip 4\p@ \@plus 2\p@ \@minus5\p@
  \belowdisplayshortskip 4\p@ \@plus 2\p@ \@minus5\p@
 \parindent 11pt %\belowdisplayskip \abovedisplayskip
  \let\@listi\@listI\baselineskip12pt
  }

\def\Hline{\noalign{\ifnum0=`}\fi{\hskip-11pt\textcolor{tabrulecolor}{\rule{\textwidth}{0.3pt}}}\futurelet\@tempa\@xhline}

\def\fulltoprule{\\[-3pt]\Hline\\[-9.5pt]}
\def\fullmidrule{\\[-9pt]\Hline\\[-9pt]}
\def\fullbotrule{\\[-8pt]\Hline\\[-10pt]}

\def\Haline{\noalign{\ifnum0=`}\fi{\hskip-11pt\textcolor{tabrulecolor}{\rule{238pt}{0.3pt}}}\futurelet\@tempa\@xhline}
\def\Hbline{\noalign{\ifnum0=`}\fi{\hskip-11pt\textcolor{tabrulecolor}{\rule{238pt}{0.25pt}}}\futurelet\@tempa\@xhline}

\def\halftoprule{\\[-3pt]\Haline\\[-8pt]}
\def\halfmidrule{\\[-9pt]\Hbline\\[-8pt]}
\def\halfbotrule{\\[-10.2pt]\Haline\\[-11.7pt]}

%\def\Hbline{\noalign{\ifnum0=`}\fi{\hskip-11pt\textcolor{tabrulecolor}{\rule{2378pt}{0.3pt}}}\futurelet\@tempa\@xhline}

\def\boxmidrule{\hspace*{-12.5pt}\rule{236.5pt}{.3pt}\\[3pt]}

\setlength\lineskip{0\p@}
\setlength\normallineskip{0\p@}


\setlength\parskip{0pt\relax}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\setcounter{topnumber}{2}
\renewcommand\topfraction{.7}
\setcounter{bottomnumber}{1}
\renewcommand\bottomfraction{.3}
\setcounter{totalnumber}{3}
\renewcommand\textfraction{.2}
\renewcommand\floatpagefraction{.5}
\setcounter{dbltopnumber}{2}
\renewcommand\dbltopfraction{.7}
\renewcommand\dblfloatpagefraction{.5}
%%%%%%%%%%%%%%just page numbers%%%%%%%%%%


\if@twoside
  \def\ps@headings{%
      \def\@oddfoot{\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lr@{}}\\[-43.4pt]
\multicolumn{2}{@{}r@{}}{\textcolor{sectioncolor}{\rule{\textwidth}{0.3pt}}}\\[4pt]
\hfill      \catchfont{\textit{Genome Biol. Evol.}}\hspace*{2.04pt} \@vol:\fpage--\lpage.\quad\@rundoi\quad\@advance
&     \catchfont{\thepage}
      \end{tabular*}}
      \def\@evenfoot{\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}l@{}}\\[-45pt]
      \multicolumn{1}{@{}c@{}}{\textcolor{sectioncolor}{\rule{\textwidth}{0.3pt}}}\\[4pt]
      \thepage\hspace*{25pt}\catchfont{\textit{Genome Biol. Evol.}}\hspace*{2.04pt} \@vol:\fpage--\lpage.\quad\@rundoi\quad\@advance\hfill
\firstfont
      \end{tabular*}}
      \def\@evenhead{\CROPMARKS\gridev\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}ll@{}}\\[3pt]
      {\lrafont\@lrh}\hfill \lrjfont \textcolor{runheadcolor}{GBE} \\[-1.5pt]
      \hskip.5pt\textcolor{sectioncolor}{\rule{\textwidth}{0.3pt}} \end{tabular*} }%
     \def\@oddhead{\CROPMARKS\gridod\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lr@{}}\\[3pt]
     {\lrafont\@rrh}\hfill  \lrjfont \textcolor{runheadcolor}{GBE}\\[-1.5pt]
     \hskip-0pt\textcolor{sectioncolor}{\rule{\textwidth}{0.3pt}} \end{tabular*}}%
        }
\else
  \def\ps@headings{%
    \let\@oddfoot\@empty
     \def\@evenhead{\CROPMARKS\gridev\thepage \centering \leftmark}%
      \def\@oddhead{\CROPMARKS\gridod\hfill{\slshape \rightmark}\hfill \thepage
      }%
     \let\@mkboth\markboth
    \def\sectionmark##1{%
      }
\fi

\def\jname#1{\let\@jn=#1
\gdef\@modeq{#1}}


\def\ps@myheadings{%
    \let\@oddfoot\@empty\let\@evenfoot\@empty
    \def\@evenhead{\thepage\hfil\slshape\leftmark}%
    \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
    \let\@mkboth\@gobbletwo
    \let\sectionmark\@gobble
    \let\subsectionmark\@gobble
    }
  \if@titlepage
  \newcommand\maketitle{\begin{titlepage}%
  \let\footnotesize\small
  \let\footnoterule\relax
  \let \footnote \thanks
  \null\vfil
  \vskip 60\p@
  \begin{center}%
    {\LARGE \@title \par}%
    \vskip 2em%
    {\large
     \lineskip 0em%.75em%
      \begin{tabular}[t]{c}%
        \@author
      \end{tabular}\par}%
      \vskip 1.5em%
    {\large \@date \par}%       % Set date in \large size.
  \end{center}\par
  \@thanks
  \vfil\null
  \end{titlepage}%
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}
\else
\newcommand\maketitle{\par
  \begingroup
   \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
     \def\@makefnmark{{{\normalfont\@thefnmark}}}%
    \long\def\@makefntext##1{\parindent 15pt\noindent %%1em to 15ptAmbika
            \hb@xt@15pt{% 12pt to 15pt Ambika
              \hss\@textsuperscript{\normalfont\@thefnmark\hskip3pt}}##1}%Ambika                   %%hskipAmbika
           \if@twocolumn
      \ifnum \col@number=\@ne
        \@maketitle
      \else
   \twocolumn[\@maketitle]%
      \fi
    \else
      \newpage
      \global\@topnum\z@   % Prevents figures from going at top of page.
      \@maketitle
    \fi
    \ifx M\@jn
        \thispagestyle{plaingbe}
        \fi
    \@thanks
  \endgroup
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
%%    \global\let\@affiliation\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}
\def\@maketitle{%
  \newpage
  \null%\firstpage
  \vspace*{-8pt}
      \noindent
  \let \footnote \thanks
    {\titlefont\textcolor{pscolor3}{\@title} \par}%
    \vskip 13pt%
    {
       {\authorfont\@author}
       }%
}
\fi

\setcounter{secnumdepth}{0}
\newcounter {part}
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]
\renewcommand\thepart {\@Roman\c@part}
\renewcommand\thesection{\@arabic\c@section}
\renewcommand\thesubsection {\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection{\thesubsection.\@arabic\c@subsubsection}
\renewcommand\theparagraph{\thesubsubsection.\@arabic\c@paragraph}
\renewcommand\thesubparagraph{\theparagraph.\@arabic\c@subparagraph}
\newcommand\part{%
   \if@noskipsec \leavevmode \fi
   \par
   \addvspace{4ex}%
   \@afterindentfalse
   \secdef\@part\@spart}

\def\@part[#1]#2{%
    \ifnum \c@secnumdepth >\m@ne
      \refstepcounter{part}%
      \addcontentsline{toc}{part}{\thepart\hspace{1em}#1}%
    \else
      \addcontentsline{toc}{part}{#1}%
    \fi
    {\parindent \z@ \raggedright
     \interlinepenalty \@M
     \normalfont
     \ifnum \c@secnumdepth >\m@ne
       \Large\bfseries \partname\nobreakspace\thepart
       \par\nobreak
     \fi
          \huge \bfseries #2%
     \markboth{}{}\par}%
    \nobreak
    \vskip 3ex
    \@afterheading}
\def\@spart#1{%
    {\parindent \z@ \raggedright
     \interlinepenalty \@M
     \normalfont
     \huge \bfseries #1\par}%
     \nobreak
     \vskip 3ex
     \@afterheading}
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Sectioning
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%


\@secpenalty = -300


\newif\if@noskipsec \@noskipsectrue


\def\@startsection#1#2#3#4#5#6{\if@noskipsec \leavevmode \fi
   \par \@tempskipa #4\relax
   \@afterindentfalse
   \ifdim \@tempskipa <\z@ \@tempskipa -\@tempskipa \@afterindentfalse\fi
   \if@nobreak \everypar{}\else
     \addpenalty{\@secpenalty}\addvspace{\@tempskipa}\fi \@ifstar
     {\@ssect{#3}{#4}{#5}{#6}}{\@dblarg{\@sect{#1}{#2}{#3}{#4}{#5}{#6}}}}

\def\@substartsection#1#2#3#4#5#6{\if@noskipsec \leavevmode \fi
   \par \@tempskipa #4\relax
   \@afterindentfalse%\vspace*{-5pt}
   \ifdim \@tempskipa <\z@ \@tempskipa -\@tempskipa \@afterindentfalse\fi
   \if@nobreak \everypar{}\else
     \addpenalty{\@secpenalty}\addvspace{\@tempskipa}\fi \@ifstar
     {\@ssect{#3}{#4}{#5}{#6}}{\@dblarg{\@subsect{#1}{#2}{#3}{#4}{#5}{#6}}}}

\def\@subsubstartsection#1#2#3#4#5#6{\if@noskipsec \leavevmode \fi
   \par \@tempskipa #4\relax
   \@afterindentfalse
   \ifdim \@tempskipa <\z@ \@tempskipa -\@tempskipa \@afterindentfalse\fi
   \if@nobreak \everypar{}\else
     \addpenalty{\@secpenalty}\addvspace{\@tempskipa}\fi \@ifstar
     {\@ssect{#3}{#4}{#5}{#6}}{\@dblarg{\@subsubsect{#1}{#2}{#3}{#4}{#5}{#6}}}}

\def\@parastartsection#1#2#3#4#5#6{\if@noskipsec \leavevmode \fi
   \par \@tempskipa #4\relax
   \@afterindentfalse
   \ifdim \@tempskipa <\z@ \@tempskipa -\@tempskipa \@afterindentfalse\fi
   \if@nobreak \everypar{}\else
     \addpenalty{\@secpenalty}\addvspace{\@tempskipa}\fi \@ifstar
     {\@ssect{#3}{#4}{#5}{#6}}{\@dblarg{\@subsubsect{#1}{#2}{#3}{#4}{#5}{#6}}}}

\def\@parasect#1#2#3#4#5#6[#7]#8{\ifnum #2>\c@secnumdepth
     \let\@svsec\@empty\else
     \refstepcounter{#1}\edef\@svsec{\csname the#1\endcsname\hskip 2em}\fi
     \@tempskipa #5\relax
      \ifdim \@tempskipa>\z@
        \begingroup #6\relax
          \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M   #8}%
        \endgroup
       \csname #1mark\endcsname{#7}\addcontentsline
         {toc}{#1}{\ifnum #2>\c@secnumdepth \else
                      \protect\numberline{\csname the#1\endcsname}\fi
                    #7}\else
        \def\@svsechd{#6\hskip #3\relax  %% \relax added 2 May 90
                   \@svsec \itshape #8\csname #1mark\endcsname
                      {#7}\addcontentsline
                           {toc}{#1}{\ifnum #2>\c@secnumdepth \else
                             \protect\numberline{\csname the#1\endcsname}\fi
                       #7}}\fi
     \@xsect{#5}}


\def\@subsubsect#1#2#3#4#5#6[#7]#8{\ifnum #2>\c@secnumdepth
     \let\@svsec\@empty\else
     \refstepcounter{#1}\edef\@svsec{\csname the#1\endcsname\hskip 4pt}\fi
     \@tempskipa #5\relax
      \ifdim \@tempskipa>\z@
        \begingroup #6\relax
          \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M
          \subsubsectionfont \textcolor{\seccolor}{#8}}%
        \endgroup
       \csname #1mark\endcsname{#7}\addcontentsline
         {toc}{#1}{\ifnum #2>\c@secnumdepth \else
                      \protect\numberline{\csname the#1\endcsname}\fi
                    #7}\else
        \def\@svsechd{#6\hskip #3\relax  %% \relax added 2 May 90
                   \@svsec \itshape #8\csname #1mark\endcsname
                      {#7}\addcontentsline
                           {toc}{#1}{\ifnum #2>\c@secnumdepth \else
                             \protect\numberline{\csname the#1\endcsname}\fi
                       #7}}\fi
     \@xsect{#5}}



\def\@subsect#1#2#3#4#5#6[#7]#8{\ifnum #2>\c@secnumdepth
     \let\@svsec\@empty\else
     \refstepcounter{#1}\edef\@svsec{\csname the #1\endcsname\hskip 1em}\fi
     \@tempskipa #5\relax
      \ifdim \@tempskipa>\z@
             \@afterheading
             \vskip 2.5pt
        \begingroup #6\relax%asdfasdf
                  \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M \subsectionfont \textcolor{\seccolor}{#8}}%
        \endgroup
       \csname #1mark\endcsname{#7}\addcontentsline
         {toc}{#1}{\ifnum #2>\c@secnumdepth \else
                      \protect\numberline{\csname the#1\endcsname}\fi
                    #7}\else
        \def\@svsechd{#6\hskip #3\relax  %% \relax added 2 May 90
                   \@svsec #8\csname #1mark\endcsname
                      {#7}\addcontentsline
                           {toc}{#1}{\ifnum #2>\c@secnumdepth \else
                             \protect\numberline{\csname the#1\endcsname}\fi
                       #7}}\fi
     \@xsect{#5}}


\def\@sect#1#2#3#4#5#6[#7]#8{\ifnum #2>\c@secnumdepth
     \let\@svsec\@empty\else
     \refstepcounter{#1}\edef\@svsec{\csname the#1\endcsname.\hskip .5em}\fi
     \@tempskipa #5\relax
      \ifdim \@tempskipa>\z@
        \begingroup #6\relax
          \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty \@M  #8\par}%
        \endgroup
       \csname #1mark\endcsname{#7}\addcontentsline
         {toc}{#1}{\ifnum #2>\c@secnumdepth \else
                      \protect\numberline{\csname the#1\endcsname}\fi
                    #7}\else
        \def\@svsechd{#6\hskip #3\relax  %% \relax added 2 May 90
                   \@svsec #8\csname #1mark\endcsname
                      {#7}\addcontentsline
                           {toc}{#1}{\ifnum #2>\c@secnumdepth \else
                             \protect\numberline{\csname the#1\endcsname}\fi
                       #7}}\fi
     \@xsect{#5}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\@xsect#1{\@tempskipa #1\relax
      \ifdim \@tempskipa>\z@
       \par \nobreak
       \vskip \@tempskipa
       \@afterheading
    \else \global\@nobreakfalse \global\@noskipsectrue
       \everypar{\if@noskipsec \global\@noskipsecfalse
                   \clubpenalty\@M \hskip -\parindent
                   \begingroup \@svsechd \endgroup \unskip
                   \hskip -#1\relax  % relax added 14 Jan 91
                  \else \clubpenalty \@clubpenalty
                    \everypar{}\fi}\fi\ignorespaces}

\def\@ssect#1#2#3#4#5{\@tempskipa #3\relax
   \ifdim \@tempskipa>\z@
     \begingroup #4\@hangfrom{\hskip #1}{\interlinepenalty \@M #5\par}\endgroup
   \else \def\@svsechd{#4\hskip #1\relax #5}\fi
    \@xsect{#3}}

\newif\if@afterindent \@afterindenttrue

\def\@afterheading{\global\@nobreaktrue
      \everypar{\if@nobreak
                   \global\@nobreakfalse
                   \clubpenalty \@M
                   \if@afterindent \else {\setbox\z@\lastbox}\fi
                 \else \clubpenalty \@clubpenalty
                    \everypar{}\fi}}


%
\def\@hangfrom#1{\setbox\@tempboxa\hbox{#1}%
      \hangindent \wd\@tempboxa\noindent\box\@tempboxa}

\def\@startsection#1#2#3#4#5#6{\if@noskipsec \leavevmode \fi
   \par \@tempskipa #4\relax
   \@afterindentfalse
   \ifdim \@tempskipa <\z@
     \@tempskipa -\@tempskipa \@afterindentfalse
   \fi
   \if@nobreak \everypar{}\else
     \addpenalty{\@secpenalty}\addvspace{\@tempskipa}\fi \@ifstar
     {\@ssect{#3}{#4}{#5}{#6}}%
     {\@dblarg{\@sect{#1}{#2}{#3}{#4}{#5}{#6}}}}
\def\@sect#1#2#3#4#5#6[#7]#8{%bbbb
     \ifnum #2>\c@secnumdepth
       \let\@svsec\@empty\else
       \refstepcounter{#1}%
       \protected@edef\@svsec{\@seccntformat{#1}}%
     \fi
     \@tempskipa #5\relax
      \ifdim \@tempskipa>\z@
        \begingroup #6\relax
          %\@hangfrom
        {\hskip -10pt\relax%
        \textcolor{\seccolor}{\@svsec}}%
                    {\interlinepenalty \@M \textcolor{\seccolor}{#8}\par}%\textcolor{\seccolor}
        \endgroup
       \csname #1mark\endcsname{#7}\addcontentsline
         {toc}{#1}{\ifnum #2>\c@secnumdepth \else
                      \protect\numberline{\csname the #1\endcsname}\fi
                   #7}\else
        \def\@svsechd{#6\hskip #3\relax  %% \relax added 2 May 90
                   \@svsec \textcolor{\seccolor}{#8}{#8}\csname #1mark\endcsname%\textcolor{\sectioncolor}
                      {#7}\addcontentsline
                           {toc}{#1}{\ifnum #2>\c@secnumdepth \else
                           \protect\numberline{\csname #1\endcsname}%
                                     \fi
                       #7}}\fi
     \@xsect{#5}}
\def\@ssect#1#2#3#4#5{\@tempskipa #3\relax
   \ifdim \@tempskipa>\z@
     \begingroup #4\@hangfrom{\hskip #1}{\interlinepenalty \@M \textcolor{\seccolor}{#5}\par}%
     \endgroup
   \else \def\@svsechd{#4\hskip #1\relax #5}\fi
    \@xsect{#3}}

\edef\seccolor{black}
\edef\sseccolor{pscolor2}
\edef\ssseccolor{pscolor3}

\newcommand\section{\@startsection{section}{1}{0pt}{10pt plus 2pt minus 2pt}{3.4pt}
    {\edef\seccolor{pscolor1}\sectionfont}}


\newcommand\subsection{\@substartsection{subsection}{2}{0pt}{10pt plus 2pt minus 2pt}{-4.4pt}
    {\edef\seccolor{pscolor2}\subsectionfont}}


\newcommand\subsubsection{\@subsubstartsection{subsubsection}{3}{0pt}{12pt plus 2pt minus 2pt}{-5pt}
    {\edef\seccolor{pscolor1}\subsubsectionfont}}

\newcommand\paragraph{\@parastartsection{paragraph}{3}{0pt}{10pt plus 2pt minus 2pt}{-5pt}
{\edef\seccolor{pscolor4}\paragraphfont}}

\newcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
           {3.25ex \@plus1ex \@minus .2ex}%
           {-1em}%
           {\normalfont\normalsize}}
\if@twocolumn
  \setlength\leftmargini  {2em}
\else
  \setlength\leftmargini  {2.5em}
\fi
\leftmargin  \leftmargini
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\if@twocolumn
  \setlength\leftmarginv  {.5em}
  \setlength\leftmarginvi {.5em}
\else
  \setlength\leftmarginv  {1em}
  \setlength\leftmarginvi {1em}
\fi
\setlength  \labelsep  {.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
\renewcommand\theenumi{\@arabic\c@enumi}
\renewcommand\theenumii{\@alph\c@enumii}
\renewcommand\theenumiii{\@roman\c@enumiii}
\renewcommand\theenumiv{\@Alph\c@enumiv}
\newcommand\labelenumi{\theenumi.}
\newcommand\labelenumii{(\theenumii)}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\newcommand\labelitemi{\textbullet}
\newcommand\labelitemii{\normalfont\bfseries \textendash}
\newcommand\labelitemiii{\textasteriskcentered}
\newcommand\labelitemiv{\textperiodcentered}
\newenvironment{qanda}
               {\list{}{\labelwidth0pt \leftmargin0pt\itemindent-\leftmargin
                        \let\makelabel\descriptionlabel
                                      }}
               {\endlist}

\newenvironment{description}
               {\list{}{\labelwidth0pt \leftmargin32pt\itemindent-\leftmargin
                        \let\makelabel\descriptionlabel
                                      }}
               {\endlist}
\newcommand*\descriptionlabel[1]{\hspace\labelsep
                                \normalfont #1\hspace*{-1pt}}
\if@titlepage
  \newenvironment{abstract}{%
      \titlepage
      \null\vfil
      \@beginparpenalty\@lowpenalty
      %\begin{center}%
       \centering \bfseries \abstractname
        \@endparpenalty\@M
      %\end{center}}%
     {\vfil\null\endtitlepage}
\else
  \newenvironment{abstract}[2]{%
      \if@twocolumn
      {}
      \else
        \small
      \fi \hskip 4pt
    \begin{minipage}{36pc} %
            \vspace*{4pt}\fontsize{9}{10}\selectfont #1\\[10pt]\keyname#2}
            {\end{minipage}\vskip 1.5em}
\fi
\if@titlepage
  \newenvironment{msc}{%
      \titlepage
      \null\vfil
      \@beginparpenalty\@lowpenalty
      %\begin{center}%
        \bfseries \mscname
        \@endparpenalty\@M
      %\end{center}}%
     {\vfil\null\endtitlepage}
\else
  \newenvironment{msc}{%
      \if@twocolumn
        \section*       {\mscname}%
      \else
        \small
                  {\hspace{-21pt}\bfseries \mscname.\vspace{\z@}}%
               \vspace{12pt}
      \fi}
      {\if@twocolumn\else\fi}
\fi
\if@titlepage
  \newenvironment{key}{%
      \titlepage
      \null\vfil
      \@beginparpenalty\@lowpenalty
      %\begin{center}%
        \itseries \keyname
        \@endparpenalty\@M
      %\end{center}}%
     {\vfil\null\endtitlepage}
\else
  \newenvironment{key}{%
      \if@twocolumn
        \section*{\keyname}%
      \else
        \small
                  {\hspace{-21pt}\bfseries \keyname.\vspace{\z@}}%
             \fi}
      {\if@twocolumn\else\fi}
\fi
\newenvironment{verse}
               {\let\\\@centercr
                \list{}{\itemsep  \z@
                        \itemindent   0em%%
                        \listparindent\itemindent
                        \rightmargin  \leftmargin
                        \advance\leftmargin 0em}%1.5em}%
                \item\relax}
               {\endlist}

\newenvironment{quotation}
               {\list{}{\listparindent 1.5em%
                        \itemindent    \listparindent
                        \rightmargin   \leftmargin
                        \parsep        \z@ \@plus\p@}%
                \item\relax}
               {\endlist}






\newbox\boxH
\newbox\boxHrunin


\newenvironment{NL}[1]{
  \setbox\boxH=\hbox{{#1}}%
 \begin{list}{\arabic{enumii}.}{%
 \usecounter{enumii}
  \topsep=5.5pt plus 0pt minus 1.5pt%
 \partopsep=0pt
 \itemsep=0pt
 \parsep=0pt
 \labelsep=4.5pt
 \labelwidth\wd\boxH\advance\labelwidth\labelsep
 \leftmargin=12pt%24pt
 \rightmargin=0pt
 \listparindent=0pt
}}{\end{list}}%





%%%%%%%%%%%%%%%%%%%%bracketNL
%
%

\newenvironment{BNL}[1]{
  \setbox\boxH=\hbox{{#1}}%
 \begin{list}{\arabic{enumii})}{%
 \usecounter{enumii}
  \topsep=6pt plus 2pt minus 1.5pt
 \partopsep=0pt
 \itemsep=2pt
 \parsep=0pt
 \labelsep=6pt
 \labelwidth\wd\boxH\advance\labelwidth\labelsep
 \leftmargin=12pt%24pt
%% \advance\leftmargin\labelsep
 \rightmargin=0pt
 \listparindent=0pt
 %\itemindent-\wd\boxH\advance\itemindent-\labelsep
% \itemindent-\labelsep
}}{\end{list}}%



%%%%%%%%%%%%%%%%%AlphabetNL
%
%

\newenvironment{ANL}[1]{
  \setbox\boxH=\hbox{{#1}}%
 \begin{list}{(\alph{enumii})}{%
 \usecounter{enumii}
  \topsep=6pt plus 2pt minus 1.5pt
 \partopsep=0pt
 \itemsep=2pt
 \parsep=0pt
 \labelsep=5pt
 \labelwidth\wd\boxH\advance\labelwidth\labelsep
 \leftmargin=28pt
 \rightmargin=0pt
 \listparindent=0pt
}}{\end{list}}%



%%%%%%%%%%%%%%roman NL
%
%

\newenvironment{RNL}[1]{
  \setbox\boxH=\hbox{{#1}}%
 \begin{list}{(\roman{enumii})}{%
 \usecounter{enumii}
  \topsep=6pt plus 2pt minus 1.5pt
 \partopsep=0pt
 \itemsep=2pt
 \parsep=0pt
 \labelsep=5pt
 \labelwidth\wd\boxH\advance\labelwidth\labelsep
 \leftmargin18pt\advance\leftmargin\wd\boxH
 \rightmargin=0pt
 \listparindent=0pt
}}{\end{list}}%



%%%%%%%%%%Double alpha NL
%
%
\newenvironment{DANL}[1]{
  \setbox\boxH=\hbox{{#1}}%
 \begin{list}{(\alph{enumii})}{%
 \usecounter{enumii}
  \topsep=6pt plus 2pt minus 1.5pt
 \partopsep=0pt
 \itemsep=2pt
 \parsep=0pt
 \labelsep=6pt
 \labelwidth\wd\boxH\advance\labelwidth\labelsep
 \leftmargin=24pt
 \rightmargin=0pt
 \listparindent=0pt
}}{\end{list}}%



\newenvironment{BL}{
  \begin{list}{$\bullet$}{\fontsize{10}{12}\selectfont%
 \topsep=6pt plus 2pt minus 1.5pt
 \partopsep=0pt
 \itemsep=6pt
 \parsep=0pt
 \labelsep=9pt
 \labelwidth\labelsep
 \leftmargin=13pt
 \rightmargin=0pt
 \listparindent=0pt}
}
{\end{list}}%



\newenvironment{UNL}{
  \begin{list}{}{
 \topsep=0pt %6pt plus 2pt minus 1.5pt
 \partopsep=0pt
 \itemsep=0pt
 \parsep=0pt
 \leftmargin=0pt
 \rightmargin=0pt
 \listparindent=0pt}
}
{\end{list}}%



\newenvironment{references}
               { %\vskip12pt
                          \fontsize{9}{11}\selectfont
        %
           \list{}
           {\vskip12pt\centering\fboxsep0pt\mbox{{\smlc References}\vspace*{3.5pt}}
            \topsep0pt
           \labelwidth\z@%
            \leftmargin 18\p@
            \itemindent-\leftmargin
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
            \fontsize{9}{11}\selectfont\parindent\z@
            \itemsep1pt
            \parskip0pt\relax%
      \sloppy\clubpenalty4000\widowpenalty4000%
      \sfcode`\.\@m%can be changed to 1000 instead of \@m for equal spacing
      %allover

          }
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist
     }





\newenvironment{quote}
               {\list{}{\rightmargin\leftmargin}%
                \item\relax}
               {\endlist}
\if@compatibility
\newenvironment{titlepage}
    {%
      \if@twocolumn
        \@restonecoltrue\onecolumn
      \else
        \@restonecolfalse\newpage
      \fi
      \thispagestyle{empty}%
      \setcounter{page}\z@
    }%
    {\if@restonecol\twocolumn \else \newpage \fi
    }
\else
\newenvironment{titlepage}
    {%
      \if@twocolumn
        \@restonecoltrue\onecolumn
      \else
        \@restonecolfalse\newpage
      \fi
      \thispagestyle{empty}%
      \setcounter{page}\@ne
    }%
    {\if@restonecol\twocolumn \else \newpage \fi
     \if@twoside\else
        \setcounter{page}\@ne
     \fi
    }
\fi
\newcommand\appendix{\par
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  \gdef\thesection{\@Alph\c@section}}
\setlength\arraycolsep{5\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
%\skip\@mpfootins = \skip\footins
\renewcommand\theequation{\arabic{equation}}
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}
\newcounter{figure}
\renewcommand \thefigure {\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\hskip 21.5pt{{\figcapheadfont{F{\scshape\fontsize{7}{7}\selectfont ig}.\hspace*{0.8pt}}}\nobreakspace\thefigure.}}%. %add
\newenvironment{figure}
               {\noindent\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
\newcounter{table}
\renewcommand\thetable{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename\nobreakspace\thetable}
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{7pt}
\setlength\belowcaptionskip{4\p@}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
\fontsize{9}{10.5}\selectfont\noindent\hspace*{-5.8pt}{\figcapheadfont\selectfont{\textcolor{figcolor}{#1}}}---{\figcapfont #2}
  \vskip\belowcaptionskip}


%%%%%%%%%%%%%%%%%%%%%%%%
%\long\def\@makecaption#1#2{%
%  \vskip\abovecaptionskip
%  \sbox  \@tempboxa{\fontsize{9}{10.5}\selectfont{\figcapheadfont\selectfont{\textcolor{figcolor}{#1}}}\figcapfont#2}
%  \ifdim \wd\@tempboxa >\hsize
%\fontsize{9}{10.5}\selectfont\noindent\hspace*{-5.8pt}{\figcapheadfont\selectfont{\textcolor{figcolor}{#1}}} %{\figcapfont #2}
%  \else
%    \global \@minipagefalse
%    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
%  \fi
%  \vskip\belowcaptionskip}

%%%%%%%%%%%%%%%%%%%%%%%%%%


\newcounter{tblcap}
\renewcommand\thetblcap{\@arabic\c@tblcap}

\newsavebox{\tabbox}
\newlength{\tablen}

\long\def\boxcaption#1#2{
\sbox{\tabbox}{#2}
\settowidth{\tablen}{\usebox{\tabbox}}
\addtolength{\tablen}{-18pt}
\refstepcounter{tblcap}
\noindent\parbox[c]{\tablen}{\fontsize{7.8}{9}\selectfont
\noindent\mbox{\textcolor{tabcolor}{\tabcapheadfont Table \thetblcap.}}
\fontsize{7.8}{9}\selectfont{\tabcapfont #1}.}\\
\noindent \fontsize{7.8}{9}\selectfont \mbox{\tbhf #2}\\[6pt]
\parbox[l]{\tablen}{\parindent 0pt\fontsize{7.5}{8.5}\tabff\selectfont}}





\long\def\tblcaption#1#2{
\sbox{\tabbox}{#2}
\settowidth{\tablen}{\usebox{\tabbox}}
\addtolength{\tablen}{-18pt}
\refstepcounter{tblcap}
\noindent\parbox[c]{\tablen}{\fontsize{7.8}{9}\selectfont
\noindent\mbox{\textcolor{tabcolor}{\tabcapheadfont Table} \textcolor{tabcolor}{\tabcapheadfont
\thetblcap}}\\[7pt]
\fontsize{7.8}{9}\selectfont{\tabcapfont #1}}\\[-.9pt]
\noindent \fontsize{7.8}{9}\selectfont \mbox{\tbhf #2}\\[6pt]
\parbox[l]{\tablen}{\parindent 0pt\fontsize{7.5}{8.5}\tabff\selectfont}}



\long\def\fulltblcaptionnotes#1#2#3{
\sbox{\tabbox}{#2}
\settowidth{\tablen}{\usebox{\tabbox}}
\addtolength{\tablen}{-18pt}
\vskip5pt\refstepcounter{tblcap}
\noindent\parbox[c]{40pc}{\fontsize{7.8}{9}\selectfont
\noindent\mbox{\textcolor{tabcolor}{\tabcapheadfont Table} \textcolor{tabcolor}{\tabcapheadfont
\thetblcap}} \\[7pt]
\fontsize{7.8}{9}\selectfont{\tabcapfont #1}}\\[-.9pt]
\noindent \fontsize{7.8}{9}\selectfont \mbox{\tbhf #2}\\[7.6pt]
\parbox[l]{\tablen}{\parindent 14pt\fontsize{7}{9}\tabff\selectfont #3}}



\long\def\tblcaptionnotes#1#2#3{
\sbox{\tabbox}{#2}
\settowidth{\tablen}{\usebox{\tabbox}}
\addtolength{\tablen}{-18pt}
\refstepcounter{tblcap}
\noindent\parbox[c]{\tablen}{\fontsize{7.8}{9}\selectfont
\noindent\mbox{\textcolor{tabcolor}{\tabcapheadfont Table} \textcolor{tabcolor}{\tabcapheadfont
\thetblcap}}\\[7pt]
\fontsize{7.8}{1}\selectfont{\tabcapfont #1}}\\[-.9pt]
\noindent \fontsize{7.8}{9}\selectfont \mbox{\tbhf #2}\\[7.6pt]
\parbox[l]{20pc}{\parindent 14pt\fontsize{7}{9}\tabff\selectfont #3}}







\long\def\firstline#1#2{\\[-8pt]\multicolumn{#2}{@{\hspace*{-6pt}}c@{}}{\rule{#1}{1pt}} \\[2.6pt]}  %[-8pt]}
\long\def\midline{\\[1.5pt]\hline\\[-.6pt]}





\def\lastline{\\[2pt]\hline\\[-10pt]}

\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}
\RequirePackage{amssymb,amsthm,graphicx,latexsym}
\usepackage{amsmath,color,stfloats,multirow,booktabs,bm}%,modnatbib}
\usepackage{cuted,frutiger,rotating}
\usepackage[switch]{tnqlineno}
%xxx sdj pifont requires the font msyr
%xxx sdj \usepackage{tnqcuted,tnqcolor,pifont}
\usepackage{tnqcuted,tnqcolor}
\usepackage{ulem,color}
\usepackage{hyperref} %xxx sdj for href
\definecolor{pscolor1}{cmyk}{.93,.56,.12,.22}
\definecolor{pscolor2}{cmyk}{.60,.80,0,.19}
\definecolor{pscolor3}{cmyk}{.60,.80,0,.19}
\definecolor{citecolor}{cmyk}{.50,.50,.0,.50}
\definecolor{pscolor4}{cmyk}{.60,.80,0,.19}
\definecolor{pscolor5}{cmyk}{.93,.56,.12,.22}
\definecolor{pscolor6}{cmyk}{0.70,1.00,0,0.4}
\normalem



\newcommand\@pnumwidth{1.55em}
\newcommand\@tocrmarg{2.55em}
\newcommand\@dotsep{4.5}
\setcounter{tocdepth}{3}
\newcommand\tableofcontents{%
    \section*{\contentsname
        \@mkboth{%
           \MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
    \@starttoc{toc}%
    }
\newcommand*\l@part[2]{%
  \ifnum \c@tocdepth >-2\relax
    \addpenalty\@secpenalty
    \addvspace{2.25em \@plus\p@}%
    \setlength\@tempdima{3em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      {\leavevmode
       \large \bfseries #1\hfil \hb@xt@\@pnumwidth{\hss #2}}\par
       \nobreak
       \if@compatibility
         \global\@nobreaktrue
         \everypar{\global\@nobreakfalse\everypar{}}%
      \fi
    \endgroup
  \fi}
\newcommand*\l@section[2]{%
  \ifnum \c@tocdepth >\z@
    \addpenalty\@secpenalty
    \addvspace{1.0em \@plus\p@}%
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
    \endgroup
  \fi}
\newcommand*\l@subsection{\@dottedtocline{2}{1.5em}{2.3em}}
\newcommand*\l@subsubsection{\@dottedtocline{3}{3.8em}{3.2em}}
\newcommand*\l@paragraph{\@dottedtocline{4}{7.0em}{4.1em}}
\newcommand*\l@subparagraph{\@dottedtocline{5}{10em}{5em}}
\newcommand\listoffigures{%
    \section*{\listfigurename}%
      \@mkboth{\MakeUppercase\listfigurename}%
              {\MakeUppercase\listfigurename}%
    \@starttoc{lof}%
    }
\newcommand*\l@figure{\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand\listoftables{%
    \section*{\listtablename}%
      \@mkboth{%
          \MakeUppercase\listtablename}%
         {\MakeUppercase\listtablename}%
    \@starttoc{lot}%
    }
\let\l@table\l@figure
\newdimen\bibindent
\setlength\bibindent{17pt}
%
%
\newenvironment{thebibliography}[1]
     %
     {\vskip23pt\vbox{\centering{\sc References}\vrule depth 9pt width 0pt
     height 0pt}
                  %
      \list{}%
              {%
           \topsep0pt%3.5
            \fontsize{7.8}{10.5}\selectfont
            \parsep0pt
            \leftmargin18pt
            \listparindent -18pt%
                         \leftmargin 18pt
                         \itemindent    \listparindent
                        \rightmargin 0pt
                        \itemsep 0pt
                          \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            %
      }}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist
     }


\newcommand\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty
\newenvironment{theindex}
               {\if@twocolumn
                  \@restonecolfalse
                \else
                  \@restonecoltrue
                \fi
                \columnseprule \z@
                \columnsep 35\p@
                \twocolumn[\section*{\indexname}]%
                \@mkboth{\MakeUppercase\indexname}%
                        {\MakeUppercase\indexname}%
                \thispagestyle{plain}\parindent\z@
                \parskip\z@ \@plus .3\p@\relax
                  \let\item\@idxitem}
               {\if@restonecol\onecolumn\else\clearpage\fi}
\newcommand\@idxitem{\par\hangindent 40\p@}
\newcommand\subitem{\@idxitem \hspace*{20\p@}}
\newcommand\subsubitem{\@idxitem \hspace*{30\p@}}
\newcommand\indexspace{\par \vskip 10\p@ \@plus5\p@ \@minus3\p@\relax}
\renewcommand\footnoterule{%
  }
\newcommand\@makefntext[1]{
    \parindent 6pt%
    \noindent
\hb@xt@6pt{\hss\@makefnmark\hskip2pt}\fontsize{8}{9.5}\selectfont #1}%
%
\setlength{\skip\footins}{15pt}%{15.9pt}%
\newcommand\contentsname{Contents}
\newcommand\listfigurename{List of Figures}
\newcommand\listtablename{List of Tables}
\newcommand\refname{\textcolor{sectioncolor}{Literature Cited}}
\newcommand\indexname{Index}
\newcommand\figurename{aa}
\newcommand\tablename{Table}
\newcommand\partname{Part}
\newcommand\appendixname{Appendix}
\newcommand\abstractname{}%Abstract}
\newcommand\mscname{\it Mathematical subject classification 2000:\ }
\def\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
%\setlength\columnsep{12pt}%{10\p@}
\setlength\columnseprule{0\p@}
\pagenumbering{arabic}
\if@twoside
\else
  \raggedbottom
\fi
\if@twocolumn
  \twocolumn
  \sloppy
  \flushbottom
\else
  \onecolumn
\fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%begin_CROPMARKS

\newdimen\markwidth
\newdimen\marklength
\newdimen\markdistance
\newdimen\trimwidth
\newdimen\trimheight
\newdimen\rhsink
\newdimen\leftmarginodd
\newdimen\leftmargineven
\newdimen\leftmarginodda
\newdimen\leftmarginevena
\newdimen\auxaaa
\newdimen\auxbbb
\newdimen\auxccc
\newdimen\auxddd
\newdimen\tempbbb
\newdimen\rhsinka

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% newcropmarksetting
\markwidth=.5pt %0.4 truept % Thickness of the mark line
\marklength=0pt%truept % Length of the mark line
\markdistance=10pt %-6 truept % Distance of the real corner to the beg. of the mark
%\trimwidth=594truept%612truept%8.5in %189mm
\trimwidth=597.5078pt
%\trimheight=792truept%792truept%11in %246mm
\trimheight=796.6771truept%786.6771%794.97truept%792truept%11in %246mm
\rhsink=27pt%20pt %58pt%5pc %%%48 truept % Sink from the top of paper to the BASE of RH
                 %  where the crops are hooked to (2pc + height of 10pt digit
\leftmarginodd=52truept%73truept%%25mm%
\leftmargineven=48pt%54truept%24mm

\leftmarginodda=52truept%73truept%25mm%34mm %%%
\leftmarginevena=50truept%24mm%29mm %%even page left margin
\rhsinka=27pt%20pt%20pt%58pt%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 \auxaaa\marklength\advance\auxaaa\markdistance
 \def\oddevenswitch{\ifodd\c@page \auxbbb=\leftmarginodd
  \else\auxbbb=\leftmargineven\fi}%

\def\oddevenswitcha{\ifodd\c@page \tempbbb=\leftmarginodda
  \else\tempbbb=\leftmarginevena\fi}%


\def\Vtrim{\VRHDW{30pt}{0pt}{.1pt}}%
\def\Htrim{\VRHDW{.1pt}{0pt}{30pt}}%

\def\ULtrim{\koo{-\markdistance}{0pt}{\zeroLR{\Htrim}}%
  \koo{0pt}{-\markdistance}{\zeroLR{\Vtrim}}}%
\def\URtrim{\koo{\markdistance}{0pt}{\zeroLL{\Htrim}}%
  \koo{0pt}{-\markdistance}{\zeroLL{\Vtrim}}}%
\def\LLtrim{\koo{-\markdistance}{0pt}{\zeroUR{\Htrim}}%
  \koo{0pt}{\markdistance}{\zeroUR{\Vtrim}}}%
\def\LRtrim{\koo{\markdistance}{0pt}{\zeroUL{\Htrim}}%
  \koo{0pt}{\markdistance}{\zeroUL{\Vtrim}}}%

\def\CROPMARKS{\oddevenswitch
  \koo{-\auxbbb}{-\rhsink}{%
     \koo{0pt}{0pt}{\ULtrim}%
     \koo{\trimwidth}{0pt}{\URtrim}%
     \koo{0pt}{\trimheight}{\LLtrim}%
     \koo{\trimwidth}{\trimheight}{\LRtrim}
     }}

\def\CROPMARKSA{\oddevenswitcha
  \koo{-\tempbbb}{-\rhsinka}{%
     \koo{0pt}{0pt}{\ULtrim}%
     \koo{\trimwidth}{0pt}{\URtrim}%
     \koo{0pt}{\trimheight}{\LLtrim}%
     \koo{\trimwidth}{\trimheight}{\LRtrim}
     }}



\def\CROPMARKS{\oddevenswitch
  \koo{-\auxbbb}{-\rhsink}{%
     \koo{0pt}{0pt}{\ULtrim}%
     \koo{\trimwidth}{0pt}{\URtrim}%
     \koo{0pt}{\trimheight}{\LLtrim}%
     \koo{\trimwidth}{\trimheight}{\LRtrim}
     }}

\def\CCROPMARKS{\oddevenswitch
 \hbox to0pt{\kern-\auxbbb\vbox to0bp{\kern-\rhsink\hbox to\trimwidth{%
 \koo{-\auxaaa}{0pt}{\VRHDW{\markwidth}{0pt}{\marklength}}%
 \koo{-\markwidth}{-\auxaaa}{\VRHDW{0pt}{\marklength}{\markwidth}}%
 \hfill
 \koo{\markdistance}{0pt}{\VRHDW{\markwidth}{0pt}{\marklength}}%
 \koo{0pt}{-\auxaaa}{\VRHDW{0pt}{\marklength}{\markwidth}}}%
 \nointerlineskip\vskip\trimheight
 \nointerlineskip\hbox to\trimwidth{%
 \koo{-\auxaaa}{0pt}{\VRHDW{0pt}{\markwidth}{\marklength}}%
 \koo{-\markwidth}{\markdistance}{\VRHDW{0pt}{\marklength}{\markwidth}}%
 \hfill
 \koo{\markdistance}{0pt}{\VRHDW{0pt}{\markwidth}{\marklength}}%
 \koo{0pt}{\markdistance}{\VRHDW{0pt}{\marklength}{\markwidth}}}\vss}\hss}}%

%\def\CROPMARKS{} %%%%%%%%% no cropmarks
\newcommand*\l@ssection[2]{%
  \ifnum \c@tocdepth >\m@ne
    \addpenalty{-\@highpenalty}%
    \vskip 6pt plus1pt minus1pt
    \begingroup
    %\def\numberline##1{\@hangfrom{\hbox to\tocsect{##1\hss}}}%
      \parindent \z@ \raggedright
        \advance\rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \relax
    \leftskip0pt
      \fontsize{11}{13}\selectfont\sffamily
    \textcolor{pscolor2}{#1\nobreak\hfill\nobreak\hb@xt@\@pnumwidth{\hss #2}}\par
      \penalty\@highpenalty
        \vskip3pt
    \endgroup
  \fi}




\def\cropmark#1{
\let\crop=#1
\ifx N\crop
\def\CROPMARKSA{}
\def\CROPMARKS{}
\fi
}


%%%%%%%%%b-KOO
\long\def\koo#1#2#3{\vbox to0pt{\hsize0pt\kern #2\hbox to0pt{\kern
#1{#3}\hss}\vss}}
\long\def\zeroCC#1{\vbox to0pt{\vss\hbox to0pt{\hss #1\hss}\vss}}
\long\def\zeroLC#1{\vbox to0pt{\vss\hbox to0pt{\hss #1\hss}}}
\long\def\zeroUC#1{\vbox to0pt{\hbox to0pt{\hss #1\hss}\vss}}
\long\def\zeroLR#1{\vbox to0pt{\vss\hbox to0pt{\hss #1}}}
\long\def\zeroCR#1{\vbox to0pt{\vss\hbox to0pt{\hss #1}\vss}}
\long\def\zeroUR#1{\vbox to0pt{\hbox to0pt{\hss #1}\vss}}
\long\def\zeroLL#1{\vbox to0pt{\vss\hbox to0pt{#1\hss}}}
\long\def\zeroCL#1{\vbox to0pt{\vss\hbox to0pt{#1\hss}\vss}}
\long\def\zeroUL#1{\vbox to0pt{\hbox to0pt{#1\hss}\vss}}
\def\VRHDW#1#2#3{\vrule height #1 depth #2 width #3}
\def\HRHDW#1#2#3{\hrule height #1 depth #2 width #3}
\def\HD#1#2{\vrule height #1pt depth #2pt width 0pt\relax}
%%%%%%%%%e-KOO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end_CROPMARKS

\let\geq=\geqslant
\let\leq=\leqslant
\let\ge=\geqslant
\let\le=\leqslant

\newenvironment{reference}
               { \vskip12pt\centering{\smlc References\vspace*{-5.5pt}}
                          \fontsize{9}{11}\selectfont\list{}{
                         \topsep-18pt
                         \listparindent -18pt%
                         \leftmargin 18pt
                         \itemindent    \listparindent
                        \rightmargin0pt
                        \itemsep 0pt
                        \parsep0pt}%
                \item\relax}
               {\endlist}
\def\struta{\vrule height 10pt width 0pt depth 0pt}
\def\strutb{\vrule height 0pt width 0pt depth 8pt}

\usepackage{abstract}
\usepackage{gbenatbib,balance,lastpage}
\usepackage{fancyhdr}
\setlength{\bibhang}{13pt}
\setlength{\bibsep}{0pt}
%\renewcommand\bibfont{9}{10}\selectfont
\renewenvironment{thebibliography}[1]{%
\vskip-2pt
 \bibsection\parindent \z@\bibpreamble\fontsize{7.8}{10.5}\selectfont\list
   {\@biblabel{\arabic{NAT@ctr}}}{\@bibsetup{#1}%
    \setcounter{NAT@ctr}{0}}%
    \ifNAT@openbib
      \renewcommand\newblock{\par}
    \else
      \renewcommand\newblock{\hskip .11em \@plus.33em \@minus.07em}%
    \fi
    \sloppy\clubpenalty4000\widowpenalty4000
    \sfcode`\.=1000\relax
    \let\citeN\cite \let\shortcite\cite
    \let\citeasnoun\cite
    \parskip0pt
 }

{\@ifundefined{NAT@sectionbib}%
    {\renewcommand\bibsection{\chapter*{\bibname
     \@mkboth{\MakeUppercase{\markleft}}{\MakeUppercase{\markright}}}}}
    {\renewcommand\bibsection{\section*{\bibname
     \ifx\@mkboth\@gobbletwo\else\markright{\MakeUppercase{\markright}}\fi}}}}
{\renewcommand\bibsection{\section*{\refname
    \@mkboth{\MakeUppercase{\markleft}}{\MakeUppercase{\markright}}}}}


\public{N}
\renewcommand{\abstracttextfont}{\fontsize{9}{10}\selectfont\rm}
\setlength{\absleftindent}{24pt}%{27.5pt}
\setlength{\absrightindent}{\absleftindent}
\jname{M}
\jno{msp106}
\gridframe{N}

\linespread{1.0}


\def\headont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{m}\fontshape{n}\fontsize{9.7}{22}\selectfont}
\def\titlefont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{bx}\fontshape{n}\fontsize{18}{23}\selectfont\raggedright}
\def\lrafont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{l}\fontshape{n}\fontsize{8.3}{12}\selectfont}
\def\firstfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{l}\fontshape{n}\fontsize{10}{12}\selectfont}
\def\lrdfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{l}\fontshape{n}\fontsize{8.3}{12}\selectfont}
\def\lrjfont{\fontencoding{OT1}\fontfamily{pmn}\fontseries{m}\fontshape{n}\fontsize{30}{30}\selectfont}
\def\authorfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{l}\fontshape{n}\fontsize{11}{13}\selectfont}
\def\afffont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{l}\fontshape{n}\fontsize{8.3}{14}\selectfont}
\def\abheadfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{eb}\fontshape{n}\fontsize{12}{14.5}\selectfont}
\def\abtextfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{l}\fontshape{n}\fontsize{9.5}{12}\selectfont}
\def\sectionfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{eb}\fontshape{n}\fontsize{12}{14.5}\selectfont}
\def\subsectionfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{bx}\fontshape{n}\fontsize{9.5}{12}\selectfont}
\def\subsubsectionfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{b}\fontshape{n}\fontsize{9.5}{12}\selectfont}
\def\paragraphfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{b}\fontshape{n}\fontsize{9.5}{12}\selectfont}
\def\copyfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{l}\fontshape{n}\fontsize{7}{9}\selectfont}
\def\catchfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{l}\fontshape{n}\fontsize{8.3}{11.2}\selectfont}
\def\tabcapheadfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{eb}\fontshape{n}\fontsize{7.8}{9}\selectfont}
\def\tabcolheadfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{bx}\fontshape{n}\fontsize{7.47}{9}\selectfont}
\def\tbhf{\fontencoding{OT1}\fontfamily{pfr}\fontseries{l}\fontshape{n}\fontsize{7.5}{10.5}\selectfont}
\def\tabbf{\fontencoding{OT1}\fontfamily{pfr}\fontseries{m}\fontshape{n}\fontsize{8.5}{10}\selectfont}
\def\tabff{\fontencoding{OT1}\fontfamily{pfr}\fontseries{l}\fontshape{n}\fontsize{6.5}{8.5}\selectfont}
\newcommand{\tabfn}[1]{\fontencoding{OT1}\fontfamily{pfr}\fontseries{l}\fontshape{n}\fontsize{5}{5}\selectfont
$^{\rm #1}$\hspace*{3.7pt}}
\def\tabcapfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{l}\fontshape{n}\fontsize{7.8}{9}\selectfont}
\def\mfigcapfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{bx}\fontshape{n}\fontsize{9}{11}\selectfont}
\def\openaccessfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{bx}\fontshape{n}\fontsize{20}{22}\selectfont}
\def\figcapheadfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{eb}\fontshape{n}\fontsize{7.8}{10.5}\selectfont}
\def\figcaptextheadfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{bx}\fontshape{n}\fontsize{7.8}{9}\selectfont}
\def\figcapfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{l}\fontshape{n}\fontsize{7.8}{9}\selectfont}
\def\keyfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{bx}\fontshape{n}\fontsize{9.5}{12}\selectfont}
\def\norfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{l}\fontshape{n}\fontsize{9.5}{12}\selectfont}

\def\dagfont{\fontencoding{OT1}\fontfamily{phv}\fontseries{m}\fontshape{n}\selectfont}

\newcommand\gbeabstract[2]{\twocolumn[
\maketitle
\begin{onecolabstract}\vspace*{-2.3pt}
\hspace*{-35pt}{\abheadfont \textcolor{sectioncolor}{Abstract}}\vspace*{13.4pt}\\
\hspace*{-25pt}
\begin{minipage}{497.9232pt}\abtextfont\noindent \color{pscolor6}{#1}\vspace*{3pt}\\
\textcolor{abscolor}{\keyfont{Key words:}} \textcolor{pscolor6}{#2}\end{minipage}
\end{onecolabstract}\vspace*{3pt}
]
}

%\newcommand\mbeheader[2]{\markboth{#1}{#2}}
\newcommand\gbestyle{%\jourinfo
\pagestyle{headings}}

\sfcode`\. 1000
\gdef\affiliation#1{\vskip-9.6pt\afffont #1}
\cropmark{Y}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%AUTHOR QUERIES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newcommand\marginlabel[1]{\mbox{}{\marginpar{\vspace*{-5.5pt}\raggedright\hspace{5pt}#1}}}

%\def\leftquery#1{\marginlabel{\hspace*{-30pt}\fontsize{8}{8}\selectfont {\color{black} {\bf {#1}}}}}

\def\leftquery#1{\marginlabel{\hspace*{30pt}\fontsize{8}{8}\selectfont {\color{black} {\bf {#1}}}}}

\def\rightquery#1{\marginpar{\hspace*{-22.25pt}\fontsize{8}{8}\selectfont {\color{black} {\bf {#1}}}}}

\def\rightrightquery#1{\marginpar{\hspace*{280pt}\fontsize{8}{8}\selectfont {\color{black} {\bf {#1}}}}}

\newcommand{\rquery}[1]{\rightquery{\hypertarget{#1}{#1}}}
\newcommand{\lquery}[1]{\leftquery{\hypertarget{#1}{#1}}}



\def\hisfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{bx}\fontshape{n}\fontsize{9.7}{11.2}\selectfont}

\def\aedfont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{bx}\fontshape{n}\fontsize{8.3}{9}\selectfont}
\def\aednamefont{\fontencoding{OT1}\fontfamily{pfr}\fontseries{l}\fontshape{n}\fontsize{8.3}{9}\selectfont}

\newcommand{\arthistory}[1]{\\[-9.5pt]{\hisfont{*Accepted:}} #1}

\newcommand{\aeditor}[1]{\vspace*{11.5pt}\noindent{\aedfont{Associate editor:}} \aednamefont{#1}}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\aadate#1{\gdef\@advance{#1}}
\def\jno#1{\gdef\@kxi{#1}}
\def\fpage{1}
\def\lpage{1}
\def\volume#1{\gdef\@vol{#1}}
\def\issue#1{\gdef\@iss{#1}}
\def\leftrunninghead#1{\gdef\@lrh{#1}}
\def\rightrunninghead#1{\gdef\@rrh{#1}}
\def\rdoi#1{\gdef\@rundoi{#1}}



\def\previewa{%
    \vbox to 0pt{\vss\noindent
    \definecolor{ticolor}{cmyk}{0,0,0,1}
{\noindent\begin{picture}(0,0)
\put(10,298){\titlefont\textbf\smash{{\hspace*{500pt}\textcolor{ticolor}{\vrule height13pt width40pt depth138pt}\hskip-31.4pt\begin{turn}{270}{\textcolor{white}{\hskip0ptResearch article}}\end{turn}}}}
\end{picture}}
}
}

\def\reviewa{%
    \vbox to 0pt{\vss\noindent
    \definecolor{ticolor}{cmyk}{.24,.90,.72,.14}
{\noindent\begin{picture}(0,0)
\put(10,298){\titlefont\textbf\smash{{\hspace*{500pt}\textcolor{ticolor}{\vrule height13pt width33pt depth140pt}\hskip-25pt\begin{turn}{270}{\textcolor{white}{\hskip0ptResearch article}}\end{turn}}}}
\end{picture}}
}
}



%\def\boxrule#1{\rule{228pt}{.5pt}#1}


%\long\def\secondpage{%
%\textwidth498pt% Use actual second page trim width
%\textheight695pt% Use actual second page trim height
%\relax
%}


\long\def\firstpage{%
\break\newpage\noindent\relax
}

\long\def\tnqpage{
\break\vfill{\pagebreak}\noindent\relax}

\def\tblhead#1{\tabcolheadfont #1}
\def\textbf#1{\fontencoding{OT1}\fontfamily{pfr}\fontseries{bx}\fontshape{n}\selectfont #1}


\newenvironment{textbox}[2]{
\setlength{\fboxrule}{1pt}
\begin{table}[t!]
\vspace*{3pt}\hspace*{-14pt}
\framebox[237pt]{\hspace*{-2.5pt}
{\begin{minipage}{211pt}\vspace*{6.3pt}
{\textbf{\textcolor{tabcolor}{#1}}}\\[-1.5pt] \boxmidrule\\
\hspace*{6pt}
#2\vspace*{11.7pt}\end{minipage}}}
\end{table}}

\endinput
%%
