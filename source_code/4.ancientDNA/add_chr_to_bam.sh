#! /usr/bin/bash

samtools view -H $1 | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - $1 > $2

echo "Finished!"

