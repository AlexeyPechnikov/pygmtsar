#!/bin/sh

text="$1"
fname="$2"

# notify user via Telegram
if [ -n "$TELEGRAM_TOKEN" ]
then
    curl \
        -F "chat_id=${TELEGRAM_CHAT_ID}" \
        -F caption="${TELEGRAM_SENDER}: ${text}" \
        -F document="@${fname}" \
        "https://api.telegram.org/bot${TELEGRAM_TOKEN}/sendDocument"
fi
