#!/bin/sh

# Exit immediately if a command exits with error
set -e

cluster=$1
num_nodes=$2

starcluster put $cluster provision_node.sh /data
starcluster sshmaster $cluster '/data/provision_node.sh'
for (( i=1; i<num_nodes; i++ )); do
    starcluster sshmaster $cluster "ssh $(printf 'node%03d' $i) \
        bash -l /data/provision_node.sh" &
done

