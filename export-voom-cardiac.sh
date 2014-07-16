rm -rf ./voom-cardiac-export &&
svn export --non-recursive svn+ssh://valinor.seas.ucla.edu/home/svn/voom/trunk voom-cardiac-export &&
svn export --non-recursive svn+ssh://valinor.seas.ucla.edu/home/svn/voom/trunk/src voom-cardiac-export/src &&
svn export --non-recursive svn+ssh://valinor.seas.ucla.edu/home/svn/voom/trunk/src/Applications voom-cardiac-export/src/Applications &&
svn export svn+ssh://valinor.seas.ucla.edu/home/svn/voom/trunk/src/Applications/ActionPotential voom-cardiac-export/src/Applications/ActionPotential &&
svn export svn+ssh://valinor.seas.ucla.edu/home/svn/voom/trunk/src/Elements voom-cardiac-export/src/Elements &&
svn export svn+ssh://valinor.seas.ucla.edu/home/svn/voom/trunk/src/Node voom-cardiac-export/src/Node &&
svn export svn+ssh://valinor.seas.ucla.edu/home/svn/voom/trunk/src/Quadrature voom-cardiac-export/src/Quadrature &&
svn export svn+ssh://valinor.seas.ucla.edu/home/svn/voom/trunk/src/Shape voom-cardiac-export/src/Shape 
svn export svn+ssh://valinor.seas.ucla.edu/home/svn/voom/trunk/configure.in.cardiac voom-cardiac-export/configure.in
svn export svn+ssh://valinor.seas.ucla.edu/home/svn/voom/trunk/src/Makefile.am.cardiac voom-cardiac-export/src/Makefile.am
svn export svn+ssh://valinor.seas.ucla.edu/home/svn/voom/trunk/src/Elements/Makefile.am.cardiac voom-cardiac-export/src/Elements/Makefile.am
