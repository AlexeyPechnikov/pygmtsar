#!/bin/sh

text="$1"
fname1="$2"
fname2="$3"

# notify user via Telegram
if [ -n "$TELEGRAM_TOKEN" ]
then
    name1=$(basename "$fname1")
    name2=$(basename "$fname2")
    json='[{"type":"document","media":"'attach://${name1}'"},{"type":"document","media":"'attach://${name2}'","caption":"${CAPTION}"}]'
    json=$(echo "$json" | sed "s|\${CAPTION}|${TELEGRAM_SENDER}: ${text}|")
    curl \
        -F "chat_id=${TELEGRAM_CHAT_ID}" \
        -F media="$json" \
        -F "${name1}=@${fname1}" \
        -F "${name2}=@${fname2}" \
        -H "Content-Type:multipart/form-data" \
        "https://api.telegram.org/bot${TELEGRAM_TOKEN}/sendMediaGroup"
fi
