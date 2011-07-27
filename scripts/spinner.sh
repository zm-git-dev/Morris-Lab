#!/bin/bash

logfile=/tmp/mylog
logsize=0
spinpause=0.10
linelen=0


# Output last line from log file.
function lastout()
{
    local line=$(tail -n 1 $logfile 2>/dev/null)
    if [[ "$line" ]]; then
        echo -n "     $line"

        # Erase any extra from last line.
        local len
        let len=$linelen-${#line}
        while [[ $len -gt 0 ]]
        do
            echo -n " "
            let len--
        done
        linelen=${#line}
    fi
}

# Output a spin character.
function spinout()
{
    local spinchar="$1"
    local sz
    local ll
    if [[ -f $logfile ]]; then
        echo -n -e "\r$spinchar"
        sleep $spinpause

        # Check for new message.
        sz=$(stat --printf '%s' $logfile 2>/dev/null)
        if [[ $sz -gt $logsize ]]; then
            lastout
            logsize=$sz
        fi
    fi
}

if [[ -f $logfile ]]; then
    logsize=$(stat --printf '%s' $logfile 2>/dev/null)
    if [[ $logsize -gt 0 ]]; then
        echo -n " "
        lastout
    fi

    while [[ -f $logfile ]]
    do
        spinout "/"
        spinout "-"
        spinout "\\"
        spinout "|"
        spinout "/"
        spinout "-"
        spinout "\\"
        spinout "|"
    done
    echo
fi
