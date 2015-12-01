% SJE: damplet.cls --- adapted from Simon Moore's style for Computer Lab.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simon's Lab Letter Class for 2001
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simon Moore, University of Cambridge, Computer Laboratory
% Please use freely (it's just a hack!)
% Source from: http://www.cl.cam.ac.uk/~swm11/letters/labletdoc/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Class Options:
%  usebwarms    - uses a black and white version of the Univ. arms
%  usecmfont    - uses computer modern font rather than Sabon (works
%                 better with xdvi
%  usetimesfont - uses times font rather than Sabon
%
%  e.g.:
%  \documentclass[usecmfont,usebwarms,usetimesfont]{lablet2001}
%
% Compulsory commands:
%  \fromwhome{your name}{qualifications}{title}{email}{phone}{fax}
%  \begin{letter}{their name}{their address}
%  \end{letter}
%
% Extra commands:
%  \lettertitle{title of letter}
%  \setdate{date string}
%  \setsigned{greeting}{signature}
%  \setsalutation{salutation}
%  \cclist{list of people\\to cc the\\message to}
%  \enclosures{list of enclosures} 
%  \reference{reference string} 
%  \bwarms              % black and white arms
%
% No name version of letter (as letter but recipents name is not added
% to the address):
%  \begin{letternoname}{their name}{their address}
%  \end{letternoname}


\ProvidesClass{damtplet}[2005/12/16]


\typeout{Damtplet: modified from Simon Moore's Computer Lab Letter Class}


\usepackage{sabon}                % Univ. font is Sabon!
\renewcommand{\ttdefault}{cmtt}   % but don't use courier for typewriter
\usepackage{color}
\usepackage{graphics}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A few new options (rest are passed through to article class).

\DeclareOption{usebwarms}{\AtBeginDocument{\bwarms}}
\DeclareOption{usetimesfont}{\typeout{Replacing Sabon with Times}
  \AtEndOfClass{\usepackage{times}}}
\DeclareOption{usecmfont}{\typeout{Replacing Sabon with Computer Modern}
  \renewcommand{\ComputerLaboratory}{\fontsize{12pt}{12pt}\bf Computer\hfill Laboratory}
  \AtEndOfClass{\renewcommand{\rmdefault}{cmr}}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

\ProcessOptions\relax

\LoadClass[11pt,a4paper]{article}



% page layout
\pagestyle{plain}
\setlength{\oddsidemargin}{0mm}
\setlength{\evensidemargin}{0mm}
    % use -31mm for official skip, -25mm for window envelopes
\setlength{\topmargin}{-25mm}
\setlength{\textwidth}{16cm}
\setlength{\textheight}{265mm}  %235mm if address in footer
%\setlength{\footskip}{25mm}

\parskip 10pt
\parindent 0.0pt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main content commands
\newcommand{\fromnamedetails}{Use \\from to speficfy your name,
position, email, phone, fax}
\newcommand{\fromquoldetails}{your qualifications?}
\newcommand{\fromjobdetails}{your position?}
\newcommand{\fromemaildetails}{email?}
\newcommand{\fromphonedetails}{phone?}
\newcommand{\fromfaxdetails}{fax?}
\newcommand{\fromwhom}[6]{
        \renewcommand{\fromnamedetails}{#1}
        \renewcommand{\fromquoldetails}{#2}
        \renewcommand{\fromjobdetails}{#3}
        \renewcommand{\fromemaildetails}{#4}
        \renewcommand{\fromphonedetails}{#5}
        \renewcommand{\fromfaxdetails}{#6}
        }

\newcommand{\addressdetails}{Use \\towhom to set the name and address}
\newcommand{\towhomname}{Name?}

\newcommand{\referencedetails}{}
\newcommand{\reference}[1]{\renewcommand{\referencedetails}{#1}}

\newcommand{\datedetails}{\today}
\newcommand{\setdate}[1]{\renewcommand{\datedetails}{#1}}

\newcommand{\salutationdetails}{Dear \towhomname,}
\newcommand{\setsalutation}[1]{\renewcommand{\salutationdetails}{#1}}

\newcommand{\yoursdetails}{Yours sincerely,}
\newcommand{\signeddetails}{\fromnamedetails}

\newcommand{\setsigned}[2]{\renewcommand{\yoursdetails}{#1}\renewcommand{\signeddetails}{#2}}

\newcommand{\encldetails}{}
\newcommand{\enclosures}[1]{\renewcommand{\encldetails}{\vspace*{10em}\textbf{Encl}:~~\begin{minipage}[t]{100mm}#1\end{minipage}}}

\newcommand{\cclistdetails}{}
\newcommand{\cclist}[1]{\renewcommand{\cclistdetails}{\vspace*{1em}cc:~~\begin{minipage}[t]{100mm}#1\end{minipage}}}

\newcommand{\cuarms}{\includegraphics{CUni3}}
\newcommand{\bwarms}{\typeout{Using B&W version of the arms}\renewcommand{\cuarms}{\includegraphics{BWUni3}}}

\newcommand{\dostamp}{}
\newcommand{\stamp}[1]{\renewcommand{\dostamp}{#1}}

%% \newcommand{\ComputerLaboratory}{\fontsize{14pt}{14pt}\bf
%% Computer\hfill Laboratory}

\newcommand{\DAMTPTitle}{\fontsize{14pt}{14pt}\bf Department of\vspace*{.7mm}\\
  Applied Mathematics and\vspace*{.7mm}\\Theoretical Physics}

\newcommand{\lettertitle}[1]{\begin{center}\bf#1\end{center}}

\newcommand{\dolettertop}{%
\setlength{\unitlength}{1mm}
\begin{picture}(160,80)(0,0)
\put(0,60){\makebox(160,20)[tl]{\parbox{\textwidth}{\bf\fromnamedetails\fromquoldetails\\\it\fromjobdetails}}}
\put(0,0){\makebox(100,35)[tl]{\begin{tabular}{@{}p{90mm}@{}}\addressdetails\end{tabular}}}
\put(0,10){\makebox(160,24.5)[tr]{\cuarms}}
%%\put(0,10){\makebox(160,24.5)[br]{\parbox{136pt}{\DAMTPTitle}}}
\put(0,1){\makebox(168,24.5)[br]{\parbox{160pt}{\DAMTPTitle}}}
\put(0,60){\makebox(160,20)[tr]{\dostamp}}
\end{picture}

\hfill\referencedetails \\
\datedetails

\vspace*{1em}

\salutationdetails
%\let\@oddfoot\swmfoot
\thispagestyle{empty}%
\enlargethispage{10mm}%
\begin{figure}[b!]\swmfoot\end{figure}
}


\renewcommand{\footnoterule}{}
\renewcommand{\thefootnote}{}

\newcommand{\swmfoot}{\fontsize{9pt}{10pt}\rm\hfill\begin{tabular}{@{}l@{~~}l@{ }l@{}}
 & & \\
\multicolumn{2}{@{}l@{}}{Centre for Mathematical Sciences}\\
\multicolumn{2}{@{}l@{}}{Wilberforce Road}\\
\multicolumn{2}{@{}l@{}}{Cambridge~~~CB3 0WA}\\
\multicolumn{2}{@{}l@{}}{U.K.}\\
 & \\
Tel: & \fromphonedetails\\
Fax: & \fromfaxdetails\\
E-mail: & \fromemaildetails\\
\end{tabular}%\let\@oddfoot\@empty
}


\newcommand{\doletterbottom}{

\begin{minipage}{150mm}%\textwidth}
\yoursdetails

\vspace{48pt}

\signeddetails

\encldetails

\cclistdetails

\end{minipage}
\cleardoublepage}


\newenvironment{letternoname}[2]
  {\renewcommand{\addressdetails}{#2}\renewcommand{\towhomname}{#1}\dolettertop}
  {\doletterbottom}

\newenvironment{letter}[2]
  {\renewcommand{\addressdetails}{#1\\#2}\renewcommand{\towhomname}{#1}\dolettertop}
  {\doletterbottom}


