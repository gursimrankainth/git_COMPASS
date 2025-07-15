#include <iostream>
#include <fstream>
#include <vector>
#include <sys/stat.h>

using namespace std;

// Function to get the last modified time of a file
time_t getFileModificationTime(const string& filePath) {
    struct stat fileInfo;
    if (stat(filePath.c_str(), &fileInfo) == 0)
        return fileInfo.st_mtime;
    return 0;
}

// Function to check if a file exists
bool fileExists(const string& filePath) {
    struct stat buffer;
    return (stat(filePath.c_str(), &buffer) == 0);
}

// Function to copy a file
void copyFile(const string& source, const string& destination) {
    ifstream src(source, ios::binary);
    ofstream dest(destination, ios::binary);
    if (src && dest) {
        dest << src.rdbuf();
        cout << "Copied: " << source << " -> " << destination << endl;
    } else {
        cerr << "Error copying: " << source << " -> " << destination << endl;
    }
}

int main() {
    // Remote server paths only
    vector<pair<string, string>> filePairs = {
        {"/afs/cern.ch/user/g/gkainth/phast/user/u970_DVCS.cc", "/afs/cern.ch/user/g/gkainth/phastPackages/git_COMPASS/u970_DVCS.cc"},
        {"/afs/cern.ch/user/g/gkainth/phast/user/UConn_Tools.h", "/afs/cern.ch/user/g/gkainth/phastPackages/git_COMPASS/UConn_Tools.h"},
        {"/afs/cern.ch/user/g/gkainth/phast/user/UConn_Tools.cc", "/afs/cern.ch/user/g/gkainth/phastPackages/git_COMPASS/UConn_Tools.cc"},
        {"/afs/cern.ch/user/g/gkainth/postPhastDVCS.py", "/afs/cern.ch/user/g/gkainth/phastPackages/git_COMPASS/postPhastDVCS.py"},
        {"/afs/cern.ch/user/g/gkainth/postPhastExclPi0.py", "/afs/cern.ch/user/g/gkainth/phastPackages/git_COMPASS/postPhastExclPi0.py"}
    };

    // Sync logic
    for (const auto& filePair : filePairs) {
        string file1 = filePair.first;
        string file2 = filePair.second;

        bool exists1 = fileExists(file1);
        bool exists2 = fileExists(file2);

        if (!exists1 && !exists2) {
            cerr << "Both files missing: " << file1 << " and " << file2 << endl;
            continue;
        }

        if (!exists1 && exists2) {
            copyFile(file2, file1);
            continue;
        }
        if (exists1 && !exists2) {
            copyFile(file1, file2);
            continue;
        }

        // Both exist
        time_t time1 = getFileModificationTime(file1);
        time_t time2 = getFileModificationTime(file2);

        if (time1 > time2) {
            copyFile(file1, file2);
        } else if (time2 > time1) {
            copyFile(file2, file1);
        } else {
            cout << "No changes: " << file1 << " and " << file2 << " are identical." << endl;
        }
    }

    return 0;
}

