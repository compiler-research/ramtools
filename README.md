ROOT scripts to convert a SAM file to a RAM (ROOT Alignment/Map) file and to work on RAM files.

  - To convert a SAM file to a RAM file do:

```bash
  $ root
  root [0] .x samtoram.C
  root [1] .q
```

 - To test read a RAM file do:

```bash
    $ root
    root [0] .x ramreader.C
    root [1] .q
```

 - To view a region, the equivalent of `samtools view bamexample.bam chr1:10150-10300`, do:

```bash
    $ root
    root [0] .x ramview.C("ramexample.root","chr1:10150-10300")
    root [1] .q
```

