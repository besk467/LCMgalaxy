/*
 GalaxyDB
 DB of galaxy file information

 Created by: Noah Oblath
 Date: Jan 8, 2013
 */

#include "GalaxyDB.hh"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

using namespace std;

GalaxyDB::GalaxyDB() :
    		galaxyFiles(),
            galaxyDataFiles(),
    		fileNotFound(""),
    		labelTokenName("Label"),
    		dirTokenName("Dir"),
            dataTokenPrefix("Data")
{
}

GalaxyDB::GalaxyDB(const string& filename, const string& ltName, const string& dirName, const string& dataPrefix) :
    		galaxyFiles(),
            galaxyDataFiles(),
    		fileNotFound(""),
    		labelTokenName(ltName),
    		dirTokenName(dirName),
            dataTokenPrefix(dataPrefix)
{
	LoadFile(filename);
}

GalaxyDB::~GalaxyDB()
{
}

bool GalaxyDB::LoadFile(const string& filename)
{
	ifstream dbFile(filename.c_str());
	if (! dbFile.is_open())
	{
		cerr << "Error: file <" << filename << "> did not open" << endl;
		return false;
	}
	string line;
	while (! dbFile.eof())
	{
		getline(dbFile, line);
		string trimmedLine = Trim(line);
		if (trimmedLine.empty()) continue;
		if (trimmedLine[0] == '#') continue;
		if (! AddGalaxy(trimmedLine))
		{
			cerr << "Error: problem adding line <" << line << "> from file <" << filename << ">" << endl;
			dbFile.close();
			return false;
		}
	}
	dbFile.close();
	return true;
}

bool GalaxyDB::DBIsReady() const
{
	return ! galaxyFiles.empty() && ! galaxyDataFiles.empty();
}

bool GalaxyDB::HasGalaxy(const std::string& galaxyLabel) const
{
	if (galaxyFiles.find(galaxyLabel) != galaxyFiles.end()) return true;
	return false;
}

bool GalaxyDB::HasGalaxyData(const std::string& galaxyDataLabel) const
{
	if (galaxyDataFiles.find(galaxyDataLabel) != galaxyDataFiles.end()) return true;
	return false;
}

const std::string& GalaxyDB::GetFile(const string& galaxyLabel, const string& fileType) const
{
	FileMap::const_iterator fmapIt = galaxyFiles.find(galaxyLabel);
	if (fmapIt == galaxyFiles.end())
	{
		return fileNotFound;
	}

	GalaxyFiles::const_iterator filesIt = fmapIt->second.find(fileType);
	if (filesIt == fmapIt->second.end())
	{
		return fileNotFound;
	}

	return filesIt->second;
}

const GalaxyDB::GalaxyFiles& GalaxyDB::GetFiles(const std::string& galaxyLabel) const
{
	FileMap::const_iterator fmapIt = galaxyFiles.find(galaxyLabel);
	if (fmapIt == galaxyFiles.end())
	{
		throw std::exception();
	}

	return fmapIt->second;
}

const std::string& GalaxyDB::GetDataFile(const string& galaxyLabel, const string& fileType) const
{
	FileMap::const_iterator fmapIt = galaxyDataFiles.find(galaxyLabel);
	if (fmapIt == galaxyFiles.end())
	{
		return fileNotFound;
	}
    
	GalaxyFiles::const_iterator filesIt = fmapIt->second.find(fileType);
	if (filesIt == fmapIt->second.end())
	{
		return fileNotFound;
	}
    
	return filesIt->second;
}

const GalaxyDB::GalaxyFiles& GalaxyDB::GetDataFiles(const std::string& galaxyLabel) const
{
	FileMap::const_iterator fmapIt = galaxyDataFiles.find(galaxyLabel);
	if (fmapIt == galaxyFiles.end())
	{
		throw std::exception();
	}
    
	return fmapIt->second;
}

const std::string& GalaxyDB::FileNotFound() const
{
	return fileNotFound;
}

bool GalaxyDB::AddGalaxy(const string& line)
{
	vector< string > tokens;
	istringstream iss(line);
	copy(istream_iterator< string >(iss),
			istream_iterator< string >(),
			back_inserter< vector<string> >(tokens));

	string name, dir;
	GalaxyFiles files;
    GalaxyFiles dataFiles;
	for (vector< string >::const_iterator tokenIt=tokens.begin(); tokenIt != tokens.end(); tokenIt++)
	{
		size_t equalPos = tokenIt->find('=');
		string tokenName = tokenIt->substr(0, equalPos);
		string tokenValue = tokenIt->substr(equalPos+1, tokenIt->length());
		if (tokenName == labelTokenName)
		{
			name = tokenValue;
		}
		else if (tokenName == dirTokenName)
		{
			dir = tokenValue;
		}
		else
		{
			string path(tokenValue);
			if (! dir.empty())
			{
				path = dir + string("/") + path;
			}
            if (tokenName.substr(0, dataTokenPrefix.size()) == dataTokenPrefix)
            {
                dataFiles.insert(GalaxyFiles::value_type(tokenName, path));
            }
            else
            {
                files.insert(GalaxyFiles::value_type(tokenName, path));
            }
		}
	}
	if (name.empty() || files.empty())
		return false;

	galaxyFiles.insert(FileMap::value_type(name, files));
    galaxyDataFiles.insert(FileMap::value_type(name, dataFiles));

	return true;
}

string GalaxyDB::Trim(const string& line, const string& whitespace)
{
	size_t begin = line.find_first_not_of(whitespace);
	if (begin == string::npos)
		return ""; // no content

	return line.substr(begin);
}

const string& GalaxyDB::LabelTokenName() const
{
	return labelTokenName;
}

void GalaxyDB::SetLabelTokenName(const string& name)
{
	labelTokenName = name;
}

const string& GalaxyDB::DirTokenName() const
{
	return dirTokenName;
}

void GalaxyDB::SetDirTokenName(const string& name)
{
	dirTokenName = name;
}

