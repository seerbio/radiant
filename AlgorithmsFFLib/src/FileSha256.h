#pragma once

#include <QByteArray>
#include <QCryptographicHash>
#include <QFile>
#include <QString>

#ifdef RADIANT_HAS_OPENSSL_SHA256
#include <openssl/evp.h>
#include <array>
#include <memory>
#endif

namespace CandidateBundleIODetail {

#ifdef RADIANT_HAS_OPENSSL_SHA256
// An unavailable digest provider selects the original Qt implementation.
// File access and read failures retain the existing empty-hash result.
inline bool tryOpenSslFileHash(const QString &path, QString *digest) {
    *digest = QString{};
    using Context = std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)>;
    Context context(EVP_MD_CTX_new(), EVP_MD_CTX_free);
    if (!context || EVP_DigestInit_ex(context.get(), EVP_sha256(), nullptr) != 1)
        return false;
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) return true;
    QByteArray buffer(1024 * 1024, Qt::Uninitialized);
    qint64 count;
    while ((count = file.read(buffer.data(), buffer.size())) > 0) {
        if (EVP_DigestUpdate(context.get(), buffer.constData(), static_cast<std::size_t>(count)) != 1)
            return false;
    }
    if (count < 0 || file.error() != QFileDevice::NoError || !file.atEnd())
        return true;
    std::array<unsigned char, EVP_MAX_MD_SIZE> bytes{};
    unsigned int size = 0;
    if (EVP_DigestFinal_ex(context.get(), bytes.data(), &size) != 1 || size != 32)
        return false;
    *digest = QString::fromLatin1(
        QByteArray(reinterpret_cast<const char *>(bytes.data()), int(size)).toHex());
    return true;
}
#endif

inline QString fileHash(const QString &path) {
#ifdef RADIANT_HAS_OPENSSL_SHA256
    QString digest;
    if (tryOpenSslFileHash(path, &digest)) return digest;
#endif
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) return {};
    QCryptographicHash hash(QCryptographicHash::Sha256);
    if (!hash.addData(&file) || file.error() != QFileDevice::NoError) return {};
    return QString::fromLatin1(hash.result().toHex());
}

} // namespace CandidateBundleIODetail
