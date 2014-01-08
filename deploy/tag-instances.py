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

    conn = boto.connect_ec2()

    reservations = conn.get_all_instances()

    active = [instance for res in conn.get_all_instances()
              for instance in res.instances
              if instance.state in set(['pending', 'running', 'stopping', 'stopped'])]

    logging.info('%d active instances', len(active))

    master_nodes = []
    slave_nodes = []
    for instance in active:
        group_names = [g.name for g in instance.groups]
        if group_names == [a.cluster_name + '-master']:
            master_nodes.append(instance)
        elif group_names == [a.cluster_name + '-slaves']:
            slave_nodes.append(instance)

    logging.info('%d master, %d slave', len(master_nodes), len(slave_nodes))

    if master_nodes:
        conn.create_tags([i.id for i in master_nodes], {'spark_node_type': 'master'})
    if slave_nodes:
        conn.create_tags([i.id for i in slave_nodes], {'spark_node_type': 'slave'})

    if slave_nodes or master_nodes:
        ids = [i.id for l in (master_nodes, slave_nodes) for i in l]
        conn.create_tags(ids, {'Owner': 'cmccoy', 'Purpose': 'b-cell-selection'})


    logging.info("Tagged nodes.")


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
