//
// Created by anichols on 11/07/2021.
//

#include "AWSStreamifier.h"

#include <QtTest/QtTest>

class AWSStreamifierTests : public QObject
{
    Q_OBJECT

public:
    AWSStreamifierTests() = default;
    ~AWSStreamifierTests() override = default;

private Q_SLOTS:

    void streamTextFileTest();
};

void AWSStreamifierTests::streamTextFileTest() {

    QSKIP("temp");


    const QString accessKeyId = "ASIA2OXTCKO6N76GY4ON";
    const QString secretAccessKey = "mAv41v/lHnowAfm2qtEshFYWr9AUPjS94sXA5PxH";
    const QString accessToken = "IQoJb3JpZ2luX2VjEDAaCXVzLXdlc3QtMiJHMEUCIAebzaY1JUx++SsyNpwXoHsu5sp05hOqWAYNBdVWKEvPAiEAjVo4ME17cGqndRU4d8WUhR8MWqWPM3G9CcOcs7/nUVQqlQMI+f//////////ARAFGgw3MTg4NDMwNDA3MDAiDPSsZzWzzlk9Gx1k6yrpAgOc5e5WRGtuq16HmxcCRIw8Mad6GFMNzdqc8ewS2oCR7w+jpbhB81qF6p+tlRUMwFDA+JqnWbvwfPNitQjbTC1duiGhpTXzVoq9AEmPy3nCsa+u+qdc2GjStqifPI/HOw2otOt1J42oR7TJCVOQ0T530Nqa9fuzb3Y/mDZIa1aH9rDWs/2pLzJbeFkDnP3sObGjzi9CE633SNEcBlywl4joX8BQS/HFKgLaRqxX6NsVm5WJdHVOi0b7d4LctxoEie1TbDZx2nn70YvaIRAmkbvNgJKYb0Ew6DrarCmJ7hGbfak+pT6oRk0laL/kBIvMC0x4BdGba0jsBQ0okozijBnWUX2SPZln+gEOyYZCF7ZEVnb7FZNM0IZOz/C8u/YafFHSJs/g6ATze2PmWhcxTQQYxRRc6Lbn7vo32pHuIglU0ml/W5bDlj2X8x2VzVDbkV9QAj9yudtCRZSru/hn8eva+o3VR6YUkI4wrbfhtAY6pgFb90LPF3ePaiBfB9sAa1sBYs1z2DdnUSVYrXCcedIh4dmZl52LSYUgv14U2btXCfbaXSYleZyMh85r+Dea2VsGQzQIZPvJGJByh1WVlzbBRLPpCQj5cwj9bWcgYQm/KPdSZQThKOc/bhpvIaRR3+02OZPafyrEANA65iFFdIb5Vsyz8cLyEzBg4OpdEpRdYgOY5SL4uIZJ4mpnxD3tyelRogVolg0p";

    AWSStreamifier awsStreamifier;
    const bool credentialsSet = awsStreamifier.setAWSCredentials(
        accessKeyId,
        secretAccessKey,
        accessToken
        );

    const QString uri = "s3://seer-experiments/EXP24001/mzml/EXP24001_2023ms1005bX1_A.raw.mzML";

    char* fileData;
    const bool streamIsValid = awsStreamifier.streamTextFile(uri, &fileData);

}

QTEST_MAIN(AWSStreamifierTests)
#include "AWSStreamifierTests.moc"
