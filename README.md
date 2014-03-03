# `startreerenaissance`

**this is alpha-quality software**

This package extends the counting renaissance methods of [Lemey et. al. 2012](http://bioinformatics.oxfordjournals.org/content/28/24/3248.short) to estimate site-specific selection from pairwise alignments of deep sequencing data under varying coverage.

Input sequences should be in [BAM format](http://samtools.sourceforge.net/SAMv1.pdf).
A site-specific HKY JSON model specification is required - these can be generated via [fit-star](https://github.com/cmccoy/fit-star).

Sequences are processed in parallel using [Apache Spark](http://spark.apache.org/).
[BEAGLE](https://code.google.com/p/beagle-lib/) is required for likelihood calculations.
See the `deploy/` directory for an example of provisioning a cluster on [EC2](http://aws.amazon.com/ec2/).

## Compiling

`startreerenaissance` is a mix of Java and Scala; built with [sbt](http://www.scala-sbt.org/).
You may build a `.jar` file with application code and dependencies using `sbt assembly`.
Unit tests are run via `sbt test`.

## License

Code in `startreerenaissance` is licensed under the GNU Public License.
Licenses of dependencies are included in `lib/`.
