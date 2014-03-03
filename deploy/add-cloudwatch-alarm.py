#!/usr/bin/env python
import argparse
import logging
import sys

import boto


def main():
    p = argparse.ArgumentParser()
    p.add_argument('cluster_name')
    p.add_argument('--dry-run', action='store_true')
    a = p.parse_args()

    ec2 = boto.connect_ec2()

    spot_instances = [instance for res in ec2.get_all_instances()
                      for instance in res.instances
                      if instance.state in set(['pending', 'running', 'stopping', 'stopped'])
                      and instance.spot_instance_request_id is not None]

    logging.info('%d spot instances', len(spot_instances))

    candidates = [i for i in spot_instances if i.tags.get('spark_cluster_name') == a.cluster_name]

    logging.info('%d candidates', len(candidates))

    cloudwatch = boto.connect_cloudwatch()

    extant_alarms = set(i.name for i in cloudwatch.describe_alarms())

    for instance in candidates:
        dimensions = {'InstanceId': instance.id}
        alarm_name = '{0}-idle-term'.format(instance.id)

        # The name of the metric to request.  This list can be retrieved by
        # calling ListMetrics
        metric_name = 'CPUUtilization'

        # The namespace of the metric.  This can also be retrieved by calling
        # ListMetrics
        actions = ['arn:aws:automate:{0}:ec2:terminate'.format(instance.region.name),
                   'arn:aws:sns:us-west-2:602821995734:cmccoy-alarm']

        metric = cloudwatch.list_metrics(dimensions=dimensions,
                                         metric_name=metric_name)
        if not metric:
            raise ValueError("Missing: " + metric_name)
        metric = metric[0]
        if alarm_name in extant_alarms:
            logging.warn("Alarm %s already exists - overwriting", alarm_name)

        # Terminate instances when average CPU < 15% for 18
        # periods of 5 minutes (an hour and a half)
        res = metric.create_alarm(name=alarm_name,
                                  comparison='<=',
                                  threshold=10,
                                  period=300,
                                  evaluation_periods=24,
                                  statistic='Average',
                                  alarm_actions=actions,
                                  unit='Percent')
        logging.info("%s - %s", alarm_name, res)

        extant_alarms.add(alarm_name)

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
