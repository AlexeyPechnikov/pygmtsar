#!/bin/sh
# For init script logs see /var/log/daemon.log
# https://console.cloud.google.com/compute/instances?project=$PROJECT_ID
# https://cloud.google.com/compute/docs/machine-types
# gcloud compute machine-types list --zones $ZONE
#  --machine-type n2-highcpu-8  n2-standard-8
#gcloud compute instances describe gmtsar
# --maintenance-policy=TERMINATE --restart-on-failure \

PROJECT_ID=""
SERVICE_ACCOUNT=""
ZONE="asia-southeast2-a"

# 32GB disk is enough for a single interferogram processing
# n2-standard-2 instance is enough for a single interferogram processing
gcloud compute instances create \
    --project "${PROJECT_ID}" \
    --tags gmtsar gmtsar \
    --boot-disk-size 32 --boot-disk-type=pd-ssd \
    --machine-type n2-standard-2 --zone="$ZONE" \
    --image-family=debian-10 --image-project=debian-cloud \
    --service-account="${SERVICE_ACCOUNT}" \
    --scopes https://www.googleapis.com/auth/cloud-platform \
    --metadata-from-file startup-script=GMTSAR.install.debian10.sh

#NAME        ZONE               MACHINE_TYPE  PREEMPTIBLE  INTERNAL_IP  EXTERNAL_IP     STATUS
# gmtsar  asia-southeast1-b  n2-highcpu-8               10.148.0.4   35.240.222.151  RUNNING

sleep 30

gcloud compute ssh \
    --project "${PROJECT_ID}" gmtsar \
    --zone "$ZONE"

# gcloud compute scp --zone "$ZONE" --compress filename gmtsar:~/srv/
