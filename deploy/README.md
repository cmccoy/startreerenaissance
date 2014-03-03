# README

Tools to provision a [Spark][spark] cluster on EC2.

Combines tweaked versions of the [Spark EC2 scripts][ec2-scripts] with the [ansible][ansible] infrastructure automation tool.

Requires python packages `ansible` and `boto`

[spark]: http://spark.incubator.apache.org/
[ec2-scripts]: http://spark.incubator.apache.org/docs/0.8.1/ec2-scripts.html
[ansible]: http://http://www.ansibleworks.com/

# TODO

Add cloud watch to terminate once idle?

    % mon-put-metric-alarm my-Alarm \ --namespace "AWS/EC2" \ -- dimensions " InstanceId=i-abc123" \--statistic Average \ --metric- name  CPUUtilization  \ --comparison-operator LessThanThreshold \
    --threshold 10  \ --period  86400 \ --evaluation-periods  4 \ -- alarm-actions
    arn:aws:automate:us-east-1:ec2:terminate
    % mon-put-metric-alarm my-Alarm \ --namespace "AWS/EC2" \ -- dimensions " InstanceId=i-abc123" \--statistic Average \ --metric- name  CPUUtilization  \ --comparison-operator GreaterThanThreshold
    \ --threshold 10  \ --period  86400 \ --evaluation-periods  4 \ --
    alarm-actions arn:aws:automate:us-east-1:ec2:terminate
