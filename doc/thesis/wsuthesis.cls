% Washington State University thesis class for LaTeX2e.
% Copyright (C) 2002, 2003 Richard Griswold

% This program may be distributed and/or modified under the conditions of the
% LaTeX Project Public License, either version 1.2 of this license or (at your
% option) any later version.  The latest version of this license is in
%   http://www.latex-project.org/lppl.txt
% and version 1.2 or later is part of all distributions of LaTeX version
% 1999/12/01 or later.
%
% This program consists of the file wsuthesis.cls

% This class file is an attempt to create a LaTeX class that conforms to the
% Washington State University dissertation and thesis guidelines of August
% 2001.  It may or MAY NOT produce documents that conform to these guidelines,
% so you must check all of your output to ensure that it does conform to the
% guidelines.  A PDF version of the guidelines are available at
% http://www.wsu.edu/~gradsch/DISSANDTHESES.pdf.
%
% This class file is in no way supported or endorsed by Washington State
% University.

% This class is based on pieces of several different classes:
%
% muthesis.cls, 19 JUL 1999, Copyright 1999 G.D. Gough and the University of
% Manchester.
% http://www.ctan.org/tex-archive/macros/latex/contrib/supported/muthesis/
%
% adfathesis.cls, 12 JUN 1994, by Stephen Harker at the Australian Defence
% Force Academy.
% http://www.ctan.org/tex-archive/macros/latex/contrib/supported/adfathesis/
%
% uwthesis.cls, 22 JUN 1995, by Jim Fox at the University of Washington.
% http://www.ctan.org/tex-archive/macros/latex/contrib/supported/uwthesis/
%
% datetime.sty, 30 OCT 2002, by Nicola Talbot.
% http://www.tug.org/tex-archive/macros/latex/contrib/other/datetime/

% What this class needs, provides, and requires
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{wsuthesis}[2002/07/30 WSU Thesis Class]
\RequirePackage{setspace}
\RequirePackage{times}
\RequirePackage{ifthen}

% Options for documentclass
% Default to ``Master of Science'' and ``Thesis'' (MSc,The)
\newcommand{\@degreetitle}{}  % Used on title page
\newcommand{\@abvtitle}{}     % Placed after name on abstract page (for Ph.D)
\newcommand{\@thesistype}{}
\DeclareOption{PhD}{
  \renewcommand{\@degreetitle}{Doctor of Philosophy}
  \renewcommand{\@abvtitle}{, Ph.D.}
}
\DeclareOption{MSc}{
  \renewcommand{\@degreetitle}{Master of Science in Computer Science}
  \renewcommand{\@abvtitle}{, M.S.}
}
\DeclareOption{The}{\renewcommand{\@thesistype}{thesis}}
\DeclareOption{Dis}{\renewcommand{\@thesistype}{dissertation}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ExecuteOptions{MSc}          % Default is Master of Science
\ExecuteOptions{The}          % Default is Thesis
\ProcessOptions

% Default to 8.5x11 paper and report style
\LoadClass[letterpaper]{report}

% Names for various section headings
\renewcommand{\appendixname}{APPENDIX}
\renewcommand{\bibname}{BIBLIOGRAPHY}
\renewcommand{\chaptername}{CHAPTER}
\renewcommand{\contentsname}{TABLE OF CONTENTS}
\renewcommand{\indexname}{INDEX}
\renewcommand{\listfigurename}{LIST OF FIGURES}
\renewcommand{\listtablename}{LIST OF TABLES}

% Disallow page breaks at hyphens (this will give some underfull vbox's,
% so an alternative is to use \brokenpenalty=100 and manually search
% for and fix such page breaks)
\brokenpenalty=10000

% Set page margins
%\hoffset        0.0in      % These values will give
%\oddsidemargin  0.0in      % 1" margins on both
%\evensidemargin 0.0in      % edges of the page with
%\textwidth      6.5in      % 8.5x11 paper.

\hoffset        0.0in      % These values will give
\oddsidemargin  0.25in      % 1.25" margins on both
\evensidemargin 0.0in      % edges of the page with
\textwidth      6.0in      % 8.5x11 paper.

\voffset        0.0in      % These values will give
\topmargin      0.0in      % a 1" top margin, about
\headheight     0.0in      % a 0.75" gap between
\headsep        0.0in      % the text and the page
\textheight     8.5in      % number, and a 0.75"
\footskip       0.75in     % bottom margin.

% Default paragraph indentation for indenting \raggedright paragraphs
\newlength{\defparindent}
\defparindent=\parindent

% Define the rest of the macros and supply default values
\newcommand{\dept}[1]{\gdef\@dept{#1}}
\newcommand{\numsig}[1]{\gdef\@numsig{#1}}
\newcommand{\chair}[1]{\gdef\@chair{#1}}
\newcommand{\submitdate}[1]{\gdef\@submitdate{#1}}
\newcommand{\copyrightyear}[1]{\gdef\@copyrightyear{#1}}
\newcommand{\@dept}{}
\newcommand{\@numsig}{3}
\newcommand{\@chair}{}
\newcommand{\@submitdate}{%
  \ifcase\the\month\or
    January\or February\or March\or April\or May\or June\or
    July\or August\or September\or October\or November\or December
  \fi
  \space \number\the\year}

% Stanford says that December theses should have the next year as the
% copyright year
\ifnum\month=12
  \@tempcnta=\year \advance\@tempcnta by 1
  \newcommand{\@copyrightyear}{\number\the\@tempcnta}
\else
  \newcommand{\@copyrightyear}{\number\the\year}
\fi

% Set these command to empty for error checking later on
\global\let\@author\@empty
\global\let\@title\@empty
\global\let\@dept\@empty
\global\let\@chair\@empty

% New ``if'' flags and default values for the flags
\newif\ifcopyright      \copyrightfalse
\newif\iftablespage     \tablespagetrue
\newif\iffigurespage    \figurespagetrue
\newif\ifcenterhdrs     \centerhdrsfalse
\newif\ifstrhdrs        \strhdrsfalse

% New booleans for formatting table of contents
\newboolean{MadePreface}   \setboolean{MadePreface}{false}
\newboolean{HaveChapter}   \setboolean{HaveChapter}{false}
\newboolean{HaveAppendix}  \setboolean{HaveAppendix}{false}

% Add ``Page'' over page column in table of contents and lists of tables
% and figures
\addcontentsline{toc}{pagehdr}{Page}
\addcontentsline{lot}{pagehdr}{Page}
\addcontentsline{lof}{pagehdr}{Page}

% This code will convert a chapter number to a string for chapter numbers from
% 1 to 99.  It based on some code by Nicola Talbot in datetime.sty.  That code
% was not used because it didn't support printing the number in all uppercase,
% and for some reason \uppercase and \MakeUppercase wouldn't work on the text.
\newcommand{\@unitstring}[1]{%
  \ifcase#1\relax
    \or ONE%
    \or TWO%
    \or THREE%
    \or FOUR%
    \or FIVE%
    \or SIX%
    \or SEVEN%
    \or EIGHT%
    \or NINE%
  \fi
}
\newcommand{\@teenstring}[1]{%
  \ifcase#1\relax
    TEN%
    \or ELEVEN%
    \or TWELVE%
    \or THIRTEEN%
    \or FOURTEEN%
    \or FIFTEEN%
    \or SIXTEEN%
    \or SEVENTEEN%
    \or EIGHTEEN%
    \or NINETEEN%
  \fi
}
\newcommand{\@tenstring}[1]{%
  \ifcase#1\relax
    \or TEN%
    \or TWENTY%
    \or THIRTY%
    \or FOURTY%
    \or FIFTY%
    \or SIXTY%
    \or SEVENTY%
    \or EIGHTY%
    \or NINETY%
  \fi
}
\newcount\@modctr
\newcommand{\@modulo}[2]{%
  \@modctr=#1\relax
  \divide \@modctr by #2\relax
  \multiply \@modctr by #2\relax
  \advance #1 by -\@modctr
}
\newcount\@strctr
\newcommand{\@numberstring}[1]{%
  \@strctr=#1\relax
  \ifnum\@strctr>99
    \PackageError{wsuthesis}{Chapter number #1 out of range}%
    {Chapter number must be less than 100}%
  \else
    \ifnum\@strctr>19
      \divide\@strctr by 10
      \@tenstring{\@strctr}%
      \@strctr=#1\relax
      \@modulo{\@strctr}{10}%
      \ifnum\@strctr>0\relax
        -\@unitstring{\@strctr}%
      \fi
    \else
      \ifnum\@strctr>9
        \@modulo{\@strctr}{10}%
        \@teenstring{\@strctr}%
      \else
        \ifnum\@strctr>0
          \@unitstring{\@strctr}%
        \else
          \PackageError{wsuthesis}{Chapter number #1 out of range}%
          {Chapter number must be greater than 0}%
        \fi
      \fi
    \fi
  \fi
}

\newcommand{\@numberchap}{%
  \@numberstring{\csname c@chapter\endcsname}
}

% Title and (optional) copyright page
% Use a single command to create both pages, since the title page must also
% carry copyright information if the document is copyrighted.
\renewcommand{\titlepage}{
  \thispagestyle{empty}
  \normalfont
  \begin{center}

    % Title in all uppercase
    \begin{spacing}{2}
    \uppercase\expandafter{\@title}
    \end{spacing}
    \vfill

    % Author in all uppercase
    By\\[\baselineskip]
    \uppercase\expandafter{\@author}
    \vfill

    % Degree title, department/school, and date
    A \expandafter{\@thesistype} submitted in partial fulfillment of\\
    the requirements for the degree of\\[\baselineskip]
    \uppercase\expandafter{\@degreetitle}\\
    \begin{spacing}{4}\end{spacing}WASHINGTON STATE UNIVERSITY\\
    \ifx\@dept\@empty
      \relax
    \else
      \expandafter{\@dept}\\[\baselineskip]
    \fi
    \uppercase\expandafter{\@submitdate}\\[\baselineskip]

    % Optional copyright and copyright page
    \ifcopyright
      \expandafter{\copyright} Copyright by \uppercase\expandafter{\@author},
        \@copyrightyear \\
      All rights reserved\newpage
      \vspace*{\fill}
      \thispagestyle{empty}
      % The copyright page is not numbered
      \addtocounter{page}{-1}
      \expandafter{\copyright} Copyright by \uppercase\expandafter{\@author},
        \@copyrightyear \\
      All rights reserved
    \fi

  \end{center}
  \newpage
}

% Signature page
\newcommand{\signaturepage}{%
  \vspace*{\stretch{1}}

  \doublespacing
  \noindent To the Faculty of Washington State University:

  % For some reason, going directly from \doublespacing to \onehalfspacing
  % gets rid of the extra double space
  \singlespacing
  \onehalfspacing

  \indent The members of the Committee appointed to examine the
  \expandafter{\@thesistype} of \uppercase\expandafter{\@author} find it
  satisfactory and recommend that it be accepted.

  \singlespacing

  \begin{flushright}
  \begin{minipage}{3in}
    \begin{center}

    \vspace{.5in}
    \hbox to 2.95in{\hfil\shortstack{\vrule width 2.9in height 0.4pt\\Thesis Committee Chair's Name (fix in wsuthesis.cls), Chair}}
    \vspace{.5in}
    \hbox to 2.95in{\hfil\shortstack{\vrule width 2.9in height 0.4pt\\Committee Member No. 1}}
    \vspace{.5in}
    \hbox to 2.95in{\hfil\shortstack{\vrule width 2.9in height 0.4pt\\Committee Member No. 2 (repeat as needed)}}

    % Counter for number of additional signature lines.
    % Cap the number at 10, since more than this will break the formatting.
    %\newcounter{sigs}
    %\setcounter{sigs}{\@numsig}
    %\ifnum \value{sigs}>10
    %  \@latex@warning{Too many signature lines: \@numsig.  Resetting to 10.}
    %  \setcounter{sigs}{10}
    %\fi
    %\ifnum \value{sigs}<1
    %  \@latex@warning{Too few signature lines: \@numsig.  Resetting to 1.}
    %  \setcounter{sigs}{1}
    %\fi
    %\addtocounter{sigs}{-1}
    %\whiledo{\value{sigs}>0}{%
    %  \vspace{.5in}
    %  \hbox to 2.95in{\hfil\shortstack{\vrule width 2.9in height 0.4pt}}
    %  \addtocounter{sigs}{-1}
    %}
    \end{center}
  \end{minipage}
  \end{flushright}

  \vspace*{\stretch{2}}
  \newpage
}

% Acknowledgement page (optional)
\newcommand{\acknowledgment}[1]{\def\@acknowledgment{#1}}
\acknowledgment{}
\global\let\@acknowledgment\@empty
\newcommand{\acknowledgmentpage}{%
  \begin{center}
  ACKNOWLEDGEMENT\\
  \end{center}

  % Want ragged right paragraphs with indentation
  { \raggedright \parindent=\defparindent \par{\@acknowledgment} }

  \addcontentsline{toc}{chapter}{ACKNOWLEDGEMENTS}
  \newpage
}

% Abstract page
\renewcommand{\abstract}[1]{\def\@abstract{#1}}
\abstract{}
\global\let\@abstract\@empty
\newcommand{\abstractpage}{%
  % Abstract page heading
  \doublespacing
  \begin{center}
  \uppercase\expandafter{\@title}\\
  \singlespacing
  Abstract\\
  \begin{spacing}{3}\end{spacing}by
    \expandafter{\@author}\expandafter{\@abvtitle}\\
  Washington State University\\
  \expandafter{\@submitdate}\\
  \end{center}

  % Name of committee chair - do not indent
  \ifx\@chair\@empty
    \relax
  \else
    \begin{spacing}{3}\noindent Chair: \expandafter{\@chair}\end{spacing}
  \fi
  \doublespacing

  % Want ragged right paragraphs with indentation
  { \raggedright \parindent=0.5in \par{\@abstract} }

  % Add this to the table of contents
  \addcontentsline{toc}{chapter}{ABSTRACT}
  \newpage
}

% Dedication page (optional)
\newcommand{\dedication}[1]{\def\@dedication{#1}}
\dedication{}
\global\let\@dedication\@empty
\newcommand{\dedicationpage}[1]{%
  \parbox[c][\textheight]{\linewidth}{
  \begin{center}
  \textbf{Dedication}\\[\baselineskip]

  \par{\@dedication}
  \end{center}
  }
  \newpage
}

% Fill in the preface
\newcommand{\preface}{%
  % We've made the preface
  \setboolean{MadePreface}{true}

  % All pages in the preface use Roman numbering
  % \raggedright
  \pagenumbering{roman}
  % \pagestyle{plain}

  % Single space title and signature page and double space remaining pages
  \singlespacing
  \titlepage
  \signaturepage
  \doublespacing

  % Make acknowledgment page, if one
  \ifx\@acknowledgment\@empty
    \relax
  \else
    \acknowledgmentpage
  \fi

  % Create abstract page and print warning if none
  \ifx\@abstract\@empty
    \relax
  \else
    \abstractpage
  \fi

  % Create table of contents and lists of tables and figures
  \tableofcontents
  \newpage
  \iftablespage
     \addvspace{10pt}
     \addcontentsline{toc}{chapter}{LIST OF TABLES}
     \listoftables
     \newpage
  \fi
  \iffigurespage
     \addvspace{10pt}
     \addcontentsline{toc}{chapter}{LIST OF FIGURES}
     \listoffigures
     \newpage
  \fi

  % Make dedication page, if one
  \ifx\@dedication\@empty
    \relax
  \else
    \dedicationpage\\
  \fi

  % Set page numbering to Arabic for the rest of the document
  \pagenumbering{arabic}
  \newpage
}

% Set chapter heading style
\renewcommand{\@makechapterhead}[1]{%
  \ifcenterhdrs
    {
      \parindent \z@ \centering \normalfont
      \ifnum \c@secnumdepth >\m@ne
        {
          \ifstrhdrs
            \textbf{\large \@chapapp \space \@numberchap }
          \else
            \textbf{\large \@chapapp \space \thechapter}
          \fi
        }
        \par\nobreak
        \vskip 8\p@
      \fi
      \interlinepenalty\@M
      \textbf{\large \uppercase\expandafter{#1}}\par\nobreak
      \vskip 10\p@
    }
  \else
    {
      \parindent \z@ \normalfont
      \ifnum \c@secnumdepth >\m@ne
        {
          \ifstrhdrs
            \textbf{\large \@numberchap \space}
          \else
            \textbf{\large \thechapter \space}
          \fi
        }
      \fi
      \textbf{\large \uppercase\expandafter{#1}}\par\nobreak
      \vskip 10\p@
    }
  \fi
}

% Set heading style for table of contents and lists of tables and figures
% Note that these headings are not bold, while chapter headings are bold
\renewcommand{\@makeschapterhead}[1]{%
  {
    \parindent \z@ \centering \normalfont
    \interlinepenalty\@M
    \large \uppercase\expandafter{#1}\par\nobreak
    \vskip 10\p@
  }
}

% Set section and paragraph style
\renewcommand{\section}{\@startsection{section}{1}{\z@}%
    {-1.5ex \@plus-1ex \@minus -.2ex}{0.8ex \@plus.2ex}%
    {\normalfont\large\raggedright}}
\renewcommand{\subsection}{\@startsection{subsection}{2}{\z@}%
    {-1.2ex \@plus -.5ex \@minus-.2ex}{0.5ex \@plus.1ex}%
    {\normalfont\normalsize\itshape\raggedright}}
\renewcommand{\subsubsection}{\@startsection{subsubsection}{3}{\z@}%
    {-1.0ex\@plus -.5ex \@minus -.2ex}{0.3ex \@plus .1ex}%
    {\normalfont\normalsize\itshape\raggedright}}
\renewcommand{\paragraph}{\@startsection{paragraph}{4}{\z@}%
    {1.0ex \@plus.5ex \@minus.2ex}{-1em}%
    {\normalfont\normalsize\itshape\raggedright}}
\renewcommand{\subparagraph}{\@startsection{subparagraph}{5}{\parindent}%
    {1.0ex \@plus.5ex \@minus .2ex}{-1em}%
    {\normalfont\normalsize\itshape\raggedright}}

% Set style for table of contents and lists of tables and figures
\newlength{\@chapwidth}%
\renewcommand{\l@chapter}[2]{%
  \ifnum \c@tocdepth >\m@ne
    \addpenalty{-\@highpenalty}%
    \vskip 1.0em \@plus\p@
    \settowidth{\@chapwidth}{\appendixname} % not \@chapapp
    \addtolength{\@chapwidth}{\@pnumwidth}
    \setlength\@tempdima{\@chapwidth}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \normalfont
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      \expandafter\uppercase{#1}\nobreak\
        \leaders\hbox{$\m@th
        \mkern \@dotsep mu\hbox{.}\mkern \@dotsep
        mu$}\hfil\nobreak\hb@xt@\@pnumwidth{\hss #2}\par
      \penalty\@highpenalty
    \endgroup
  \fi
}

% Define a new table of contents ``comment'' line so that lines such as
% ``CHAPTER'' and ``APPENDIX'' can be added without page numbers.  A comment
% is at the same level as a chapter.  Use '\flushleft \parindent' so the
% comment will be alighned properly.
\newcommand{\l@comment}[2]{%
  \ifnum \c@tocdepth >\m@ne
    \vskip 1.0em \flushleft \parindent \z@ #1 \par
  \fi
}

% Code to add ``Page'' over page number column in table of contents and lists
% of tables and figures.  Use \par at the end so that the first entry on the
% lists of tables and figures will be aligned correctly.  Don't use \flushleft,
% or the entry won't be aligned correctly.
\newcommand{\l@pagehdr}[2]{%
  \ifnum \c@tocdepth >\m@ne
    \vskip 1.0em \flushright #1 \par
  \fi
}

\def\@chapter[#1]#2{%
  % Code to add ``CHAPTER'' to table of contents before listing any chapters.
  \ifthenelse{\equal\@chapapp\chaptername}{
    \ifthenelse{\boolean{HaveChapter}}{}{
      \addcontentsline{toc}{comment}{CHAPTER}
      \setboolean{HaveChapter}{true}
    }
  }{}

  % Code to add ``APPENDIX'' to table of contents before listing appendicies,
  % and to add the appendix half-title page.
  \ifthenelse{\equal\@chapapp\appendixname}{
    \ifthenelse{\boolean{HaveAppendix}}{}{
      \addcontentsline{toc}{comment}{APPENDIX}

      % Insert appendix half-title page
      \vspace*{2in}
      \thispagestyle{empty}
      \begin{center}
      \large\textbf{APPENDIX}
      \end{center}
      \vspace*{\fill}
      \newpage

      \setboolean{HaveAppendix}{true}
    }
  }{}

  % Code to add the table of contents listing
  \ifnum \c@secnumdepth >\m@ne
    \refstepcounter{chapter}%
    % Uncomment this to print the chapter number to stdout when LaTeX runs
    %\typeout{\@chapapp\space\thechapter.}%
    \addcontentsline{toc}{chapter}%
    {\parindent=\defparindent \indent\thechapter.\space#1}%
  \else
    \addcontentsline{toc}{chapter}{#1}%
  \fi
  \chaptermark{#1}%
  \addtocontents{lof}{\protect\addvspace{10\p@}}%
  \addtocontents{lot}{\protect\addvspace{10\p@}}%
  \@makechapterhead{#2}%
  \@afterheading
}

% Modify bibliography to add entry in table of contents
\def\textofBIB#1{\gdef\@textofBIB{#1}}  \textofBIB{BIBLIOGRAPHY}
\def\thebibliography#1{\chapter*{\textbf{\@textofBIB}\@mkboth{}{}}
  \addcontentsline{toc}{chapter}{\protect\@textofBIB}
  \list
  {[\arabic{enumi}]}{\settowidth\labelwidth{[#1]}\leftmargin\labelwidth
  \advance\leftmargin\labelsep
  \usecounter{enumi}}
  \def\newblock{\hskip .11em plus .33em minus .07em}
  \sloppy\clubpenalty4000\widowpenalty4000
  \sfcode`\.=1000\relax
}

% Issue all of the warnings at the end of the document where they'll be seen
\AtEndDocument{%
  % Issue warnings for missing required information
  \ifthenelse{\boolean{MadePreface}}{}{
    \@latex@warning
      {-=-=- No preface! Use the \noexpand\preface command -=-=-}
  }
  \ifx\@title\@empty
    \@latex@warning
      {-=-=- No title! Use the \noexpand\title command -=-=- }
  \fi
  \ifx\@author\@empty
    \@latex@warning
      {-=-=- No author! Use the \noexpand\author command -=-=- }
  \fi
  \ifx\@abstract\@empty
    \@latex@warning
      {-=-=- No abstract! Use the \noexpand\abstract command -=-=- }
  \fi
  \ifx\@dept\@empty
    \@latex@warning
      {-=-=- No department! Use the \noexpand\dept command -=-=- }
  \fi
  \ifx\@chair\@empty
    \@latex@warning
      {-=-=- No chair! Use the \noexpand\chair command -=-=- }
  \fi
}
