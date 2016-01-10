#!/bin/bash

rm -rf kermit.ply

../bin/emvs kermit/pmvs/ kermit.ply | tee kermit-emvs.log

