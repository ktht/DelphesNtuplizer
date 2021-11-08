#!/bin/bash

# source this file

BIN_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"/bin

case :$PATH: in
  *:$BIN_DIR:*) export PATH ;;  # make sure to export
  *) export PATH=$BIN_DIR:$PATH ;;
esac

