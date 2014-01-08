#!/bin/sh

set -e
set -u


SPARK_HOME=/root/spark
SPARK_JAR=$SPARK_HOME/assembly/target/scala-2.9.3/spark-assembly_2.9.3-0.8.1-incubating-hadoop1.0.4.jar
STAR_JAR=startreerenaissance-assembly-0.1.jar
MAIN=org.fhcrc.matsen.startree.spark.StarTreeSpark
MASTER=spark://$(GET http://169.254.169.254/latest/meta-data/public-hostname):7077

for f in *.bam; do
  java -Djava.library.path=/usr/local/lib \
    -cp $SPARK_JAR:$STAR_JAR \
    $MAIN \
    $MASTER \
    --jar-path $STAR_JAR \
    --executor-memory 4g \
    --sample \
    --parallelism 128 \
    --bucket startreerenaissance \
    --prefix $(basename $f .bam) \
     04-A-M_IGHV_model.json ighvdj.fasta $f
   break
done
