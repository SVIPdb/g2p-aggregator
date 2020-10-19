#!/usr/bin/env bash

ONCOKB_CREDS=".oncokb_apikey"

if [ ! -f ${ONCOKB_CREDS} ]; then
    echo "OncoKB requires an account to access their data; first, register at https://www.oncokb.org/account/register"
    echo "Once you have registered, you can view your API key at https://www.oncokb.org/account/settings"

    if [[ -t 0 ]]; then
        read -p 'API Key: ' ONCOKB_APIKEY
        echo
        echo "${ONCOKB_APIKEY}" > ${ONCOKB_CREDS}
        echo "Your input has been stored to ${ONCOKB_CREDS}, which you should treat as sensitive."
    else
        echo "Run the following command to store your API key and re-run this script:"
        echo "echo '<API Key>' > ${ONCOKB_CREDS}"
        echo "(This file will contain your credentials, so treat it as sensitive.)"
        exit 1
    fi
fi
