#!/bin/bash


set -e

echo "Running test"
cargo run --bin dietmc -- --text --relation zkfltoolbox/dct/dct.rel --instance zkfltoolbox/dct/dct.type0.ins --connection-addr 127.0.0.1:7876 --witness zkfltoolbox/dct/dct.type0.wit &
pid_prv=$!

cargo run --bin dietmc -- --text --relation zkfltoolbox/dct/dct.rel --instance zkfltoolbox/dct/dct.type0.ins --connection-addr 127.0.0.1:7876 &
pid_vrf=$!

wait $pid_prv
wait $pid_vrf

echo "All end-to-end tests passed, yay!"