# README

Tools to provision a [Spark][spark] cluster on EC2.

Combines tweaked versions of the [Spark EC2 scripts][ec2-scripts] with [Ansible][ansible] for infrastructure automation.

Requires python packages `ansible` and `boto`. Install via:

    pip install ansible boto

[spark]: http://spark.incubator.apache.org/
[ec2-scripts]: http://spark.incubator.apache.org/docs/0.8.1/ec2-scripts.html
[ansible]: http://http://www.ansibleworks.com/
