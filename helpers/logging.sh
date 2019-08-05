#!/usr/bin/env bash
#
# Log related functions will go here

# Log related options
LOG="script_log.txt"        # Name of the log file
LOG_ENABLED=0               # 1=Yes and 0=NO (Add to Log File)
LOG_DISPLAY=1               # 1=Yes and 0=No (Display on Terminal)
RECORD="record_details.txt" # Name of the record file
# Log Function

log() {
  # Log the details of it is enabled
  if [[ ${LOG_ENABLED} -eq 1 ]]; then
    echo "RLOG $(date +%y-%m-%d" "%H:%M:%S) :$1" >>"${LOG}"
  fi
  # Display if display is enabled
  if [[ ${LOG_DISPLAY} -eq 1 ]]; then
    echo "RLOG $(date +%y-%m-%d" "%H:%M:%S) :$1"
  fi
}

record() {
  echo "RLOG $(date +%y-%m-%d" "%H:%M:%S) :$1" >>"${RECORD}"
}
