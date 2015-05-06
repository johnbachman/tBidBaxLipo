#!/bin/sh

# Exit immediately if a command exits with error
set -e

cluster=$1
num_nodes=$2

starcluster put $cluster provision_node.sh /home/sgeadmin
starcluster sshmaster $cluster '/home/sgeadmin/provision_node.sh'
for (( i=1; i<num_nodes; i++ )); do
    starcluster sshnode $cluster $(printf 'node%03d' $i) \
        '/home/sgeadmin/provision_node.sh' &
done

