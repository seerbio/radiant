//
// Created by Drucifer on 11/24/2021.
//

#ifndef COMMANDLINEPARSER_H
#define COMMANDLINEPARSER_H

#include <QCommandLineParser>

// TODO make a base class out of this and subclass this file in each app.
class CommandLineParser  : public QCommandLineParser
{
    Q_DECLARE_TR_FUNCTIONS(CommandLineParser)

public:

    struct CliParameters{
        QString mzMLFilePath;
    };

    CommandLineParser();
    ~CommandLineParser() = default;

    bool validateArguments(const QStringList &args);

    [[nodiscard]] const CliParameters &getCliParams() const;


private:

    CliParameters m_cliParams;

};

#endif //COMMANDLINEPARSER_H
