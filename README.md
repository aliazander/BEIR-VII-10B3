# BEIR-VII-10B3
Steps to recreate BEIR VII figure 10B-3

The .csv file includes "Life shortening data from Tables 1, 2, and 3 of Storer and others (1979)" as stated in BEIR VII figure 10B-3 figure caption. I used a file created by Benjamin Haley. After analyzing the data and comparing it to figure 10B-3, it was clear that not all of the data was used. Information not listed in BEIR VII that we believe to be true of BEIR VII's analysis:

- Only use female mice
- Only use RFM mice
- If doses are equal (even if dose rate is not equal) average values
- Linear regression is weighted by 1/(standard deviation)^2

