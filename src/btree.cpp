/***
 * @author See Contributors.txt for code contributors and overview of BadgerDB.
 *
 * @section LICENSE
 * Copyright (c) 2012 Database Group, Computer Sciences Department, University of Wisconsin-Madison.
 */

#include "btree.h"
#include "filescan.h"
#include "exceptions/bad_index_info_exception.h"
#include "exceptions/bad_opcodes_exception.h"
#include "exceptions/bad_scanrange_exception.h"
#include "exceptions/no_such_key_found_exception.h"
#include "exceptions/scan_not_initialized_exception.h"
#include "exceptions/index_scan_completed_exception.h"
#include "exceptions/file_not_found_exception.h"
#include "exceptions/end_of_file_exception.h"


//#define DEBUG

namespace badgerdb
{

// -----------------------------------------------------------------------------
// BTreeIndex::BTreeIndex -- Constructor
// -----------------------------------------------------------------------------

BTreeIndex::BTreeIndex(const std::string & relationName,
		std::string & outIndexName,
		BufMgr *bufMgrIn,
		const int attrByteOffset,
		const Datatype attrType)
{
    // Add your code below. Please do not remove this line.
	bufMgr = bufMgrIn;
	//this.attrByteOffset -> attrByteOffset;
	attributeType = attrType;

	std::ostringstream idxStr;
	idxStr << relationName << "." << attrByteOffset ;
	std::string indexName = idxStr.str();

	Page *headerPage, *rootPage, *leafPage;
	IndexMetaInfo* meta;
	try {
		file = new BlobFile(indexName, false); // Checks if the file exists

		headerPageNum = file -> getFirstPageNo();
		bufMgr -> readPage(file, headerPageNum, headerPage); // Reads the first page (meta)
		rootPageNum = ((IndexMetaInfo*) headerPage) -> rootPageNo; // Held in the meta page

		bufMgr -> unPinPage(file, headerPageNum, false); 
		bufMgr -> unPinPage(file, rootPageNum, false);
	}
	catch (FileNotFoundException e) {
		bufMgr -> allocPage(file, headerPageNum, headerPage); // Allocate header page
    	bufMgr -> allocPage(file, rootPageNum, rootPage); // Allocate root page

		meta = (IndexMetaInfo*) headerPage; // Fills in meta info
		strcpy(meta -> relationName, relationName.c_str());
		meta -> attrByteOffset = attrByteOffset;
		meta -> attrType = attrType;
		meta -> rootPageNo = rootPageNum;
		
		LeafNodeInt *root = (LeafNodeInt *)rootPage; // Initializes the root
    	root -> rightSibPageNo = 0;

		bufMgr -> unPinPage(file, headerPageNum, false); 
		bufMgr -> unPinPage(file, rootPageNum, false); 

		// Scans Records
		FileScan scan = FileScan(relationName, bufMgr);
		scanExecuting = true;
		RecordId rid;
		try {
			while(1) {
				scan.scanNext(rid);
				std::string fileRecord = scan.getRecord();
				const char *record = fileRecord.c_str();
				insertEntry(record + attrByteOffset, rid);
			}
		}
		catch (EndOfFileException e) {
			bufMgr -> flushFile(file);
		}
	}
}


// -----------------------------------------------------------------------------
// BTreeIndex::~BTreeIndex -- destructor
// -----------------------------------------------------------------------------

BTreeIndex::~BTreeIndex()
{
    // Add your code below. Please do not remove this line.
	scanExecuting = false;
	bufMgr -> flushFile(file);
	delete file;
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

 /**
	 * Insert a new entry using the pair <value,rid>. 
	 * Start from root to recursively find out the leaf to insert the entry in. The insertion may cause splitting of leaf node.
	 * This splitting will require addition of new leaf page number entry into the parent non-leaf, which may in-turn get split.
	 * This may continue all the way upto the root causing the root to get split. If root gets split, metapage needs to be changed accordingly.
	 * Make sure to unpin pages as soon as you can.
   * @param key			Key to insert, pointer to integer/double/char string
   * @param rid			Record ID of a record whose entry is getting inserted into the index.
	**/
void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
{
    // Add your code below. Please do not remove this line.

	// Create pair
	RIDKeyPair<const void*> newEntry;
	newEntry.set(rid, key);

}

void BTreeIndex::insertRoot(){
	
}

void BTreeIndex::insertLeaf(LeafNodeInt *leafNode, RIDKeyPair<int> ridKey) // Should use parameters like this
{
	if (leafNode -> ridArray[0].page_number == 0) { // Page is empty
		leafNode -> keyArray[0] = ridKey.key;
		leafNode -> ridArray[0] = ridKey.rid;    
	}
	else {
		int i = INTARRAYLEAFSIZE; // Number of keys in the leaf node
		while((leafNode -> ridArray[i-1].page_number == 0) && i > 0) { // Gets to the end of leafNode
			i--;
		}

		while((leafNode -> keyArray[i-1] > ridKey.key) && i > 0) { // Shifts the previous ridKey pairs
			leafNode -> keyArray[i] = leafNode -> keyArray[i-1];
			leafNode -> ridArray[i] = leafNode -> ridArray[i-1];
			i--;
		}

		leafNode -> keyArray[i] = ridKey.key; // Inserts ridKey to the leafNode
		leafNode -> ridArray[i] = ridKey.rid;
	}
}

void BTreeIndex::insertNonLeaf(NonLeafNodeInt *node, PageKeyPair<int> keyPage) {
	int i = INTARRAYNONLEAFSIZE;
	while((node -> pageNoArray[i] == 0) && i > 0) {
		i--;
	}
	while((node -> keyArray[i-1] > keyPage.key) && i > 0) {
		node -> pageNoArray[i+1] = node -> pageNoArray[i]; // TODO not 100% sure if this is the right way to shift pageNo
		node -> keyArray[i] = node -> pageNoArray[i-1]; // Shifts the previous keyPage pairs
		i--;
	}

	node -> pageNoArray[i+1] = keyPage.pageNo; // Inserts keyPage pair to the node
	node -> keyArray[i] = keyPage.key;
}

// -----------------------------------------------------------------------------
// BTreeIndex::startScan
// -----------------------------------------------------------------------------

void BTreeIndex::startScan(const void* lowValParm,
				   const Operator lowOpParm,
				   const void* highValParm,
				   const Operator highOpParm)
{
    // Add your code below. Please do not remove this line.
}

// -----------------------------------------------------------------------------
// BTreeIndex::scanNext
// -----------------------------------------------------------------------------

void BTreeIndex::scanNext(RecordId& outRid) 
{
    // Add your code below. Please do not remove this line.
}

// -----------------------------------------------------------------------------
// BTreeIndex::endScan
// -----------------------------------------------------------------------------
//
void BTreeIndex::endScan() 
{
    // Add your code below. Please do not remove this line.
}

}
