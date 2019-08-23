Protein coordinate interpolation in torsional space

This is a VERY old piece of code that I once wrote to perform interpolation
of very sparsely saved trajectories of proteins, but doing it in torsional
space to create output trajectories that are both smooth and realistic.

There are some horrible kludge hacks (or, as we would say today: "AI")
to handle the cases where torsions flip close to 180 degrees in order to
avoid the entire chain rotating 360 degrees.

It is not only possible, but virtually guaranteed that the code has some
horrible bugs that might blow up in your face, but since it's open source
you can fix those ;-)

Some sunny day we might want to create a nice Gromacs tool to do this.
