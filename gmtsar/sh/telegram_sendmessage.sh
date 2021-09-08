#!/bin/sh

text="$1"

# notify user via Telegram
if [ -n "$TELEGRAM_TOKEN" ]
then
    curl \
        -F "chat_id=${TELEGRAM_CHAT_ID}" \
        -F text="${TELEGRAM_SENDER}: ${text}" \
        "https://api.telegram.org/bot${TELEGRAM_TOKEN}/sendMessage"
fi
