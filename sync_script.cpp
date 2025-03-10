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
    return 0; // Return 0 if file does not exist
}

// Function to copy a file
void copyFile(const string& source, const string& destination) {
    ifstream src(source, ios::binary);
    ofstream dest(destination, ios::binary);
    if (src && dest) {
        dest << src.rdbuf();
        cout << "Updated: " << destination << " from " << source << endl;
    } else {
        cerr << "Error copying: " << source << " -> " << destination << endl;
    }
}

int main() {
    char mode;
    cout << "Are you running this locally or remotely? (l/r): ";
    cin >> mode;

    vector<pair<string, string>> filePairs;

    if (mode == 'l') {
        // Local paths
        filePairs = {
            {"/Users/gursimran/cern/u97_DVCS.cc", "/Users/gursimran/cern/phastPackages/git_COMPASS/u97_DVCS.cc"},
            {"/Users/gursimran/cern/phast.8.032/user/u970_DVCS.cc", "/Users/gursimran/cern/phastPackages/git_COMPASS/u970_DVCS.cc"},
            {"/Users/gursimran/cern/phast.8.032/user/UConn_Tools.h", "/Users/gursimran/cern/phastPackages/git_COMPASS/UConn_Tools.h"}
        };
    } else if (mode == 'r') {
        // Remote server paths
        filePairs = {
            {"/afs/cern.ch/user/g/gkainth/phast/user/u97_DVCS.cc", "/afs/cern.ch/user/g/gkainth/phastPackages/git_COMPASS/u97_DVCS.cc"},
            {"/afs/cern.ch/user/g/gkainth/phast/user/u970_DVCS.cc", "/afs/cern.ch/user/g/gkainth/phastPackages/git_COMPASS/u970_DVCS.cc"},
            {"/afs/cern.ch/user/g/gkainth/phast/user/UConn_Tools.h", "/afs/cern.ch/user/g/gkainth/phastPackages/git_COMPASS/UConn_Tools.h"}
        };
    } else {
        cerr << "Invalid input. Use 'l' for local or 'r' for remote." << endl;
        return 1;
    }

    // Iterate over file pairs and sync them based on modification time
    for (const auto& filePair : filePairs) {
        string file1 = filePair.first;
        string file2 = filePair.second;

        time_t time1 = getFileModificationTime(file1);
        time_t time2 = getFileModificationTime(file2);

        if (time1 > time2) {
            copyFile(file1, file2);  // Editable script is newer
        } else if (time2 > time1) {
            copyFile(file2, file1);  // Backup is newer
        } else {
            cout << "No changes: " << file1 << " and " << file2 << " are identical." << endl;
        }
    }

    return 0;
}

