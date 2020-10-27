#!/usr/bin/env bash

COSMIC_CREDS=".cosmic_creds"

if [ ! -f ${COSMIC_CREDS} ]; then
    echo "COSMIC datafiles require an account to download; first, register at https://cancer.sanger.ac.uk/cosmic/register"

    if [[ -t 0 ]]; then 
        echo "Once registered, enter your credentials below so they can be used to retrieve your file"
        read -p 'Email: ' COSMIC_USERNAME
        read -sp 'Password: ' COSMIC_PASSWORD
        echo
        echo "${COSMIC_USERNAME}:${COSMIC_PASSWORD}" | base64 > ${COSMIC_CREDS}
        echo "Your input has been stored to ${COSMIC_CREDS}, which you should treat as sensitive."
    else
        echo "Once registered, run the following command to store your credentials and re-run this script:"
        echo "echo 'email@example.com:mycosmicpassword' | base64 > ${COSMIC_CREDS}"
        echo "(This file will contain your credentials, so treat it as sensitive.)"
        exit 1
    fi
fi
