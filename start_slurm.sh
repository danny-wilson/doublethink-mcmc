#!/bin/bash
# Ensure IPv4 localhost is in /etc/hosts
grep -q "^127.0.0.1.*localhost" /etc/hosts || echo "127.0.0.1 localhost" >> /etc/hosts

# Ensure Munge socket dir exists with correct permissions
mkdir -p /run/munge
chown munge:munge /run/munge
chmod 0711 /run/munge

# Start daemons
/etc/init.d/munge start
slurmctld &
slurmd &

# Drop into interactive shell
exec bash