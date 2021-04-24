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

void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
{
    // Add your code below. Please do not remove this line.
}

void BTreeIndex::insertLeaf(LeafNodeInt *leafNode, RIDKeyPair<int> ridKey) // Should use parameters like this
{
	if (leafNode -> ridArray[0].page_number == 0) { // Page is empty
		leafNode -> keyArray[0] = ridKey.key;
		leafNode -> ridArray[0] = ridKey.rid;    
	}
	else {
		int i = INTARRAYLEAFSIZE; // Number of keys in the leaf node // Maybe change to leafOccupancy ?
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

void BTreeIndex::insertNonLeaf(NonLeafNodeInt *node, PageKeyPair<int> keyPage)
{
	int i = INTARRAYNONLEAFSIZE; // Maybe change to nodeOccupancy ?
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

  PageKeyPair<int> BTreeIndex::splitLeafNode(LeafNodeInt *oldLeafNode, PageId currPageId, RIDKeyPair<int> ridPair)
  {
	// allocate new leaf sibling
	Page* sibling;
	PageId siblingId; // new leaf sibling page number
	bufMgr -> allocPage(file, siblingId, sibling);
	// convert to proper structure
	LeafNodeInt* siblingNode = (LeafNodeInt*) sibling;

	// current leaf node gets right sibling page number added to it
	if (oldLeafNode -> rightSibPageNo != 0) {
		siblingNode -> rightSibPageNo = oldLeafNode -> rightSibPageNo;
	}
	oldLeafNode -> rightSibPageNo = siblingId;

	// split current leaf into two seperate leaves
	for (int i = 0; i < INTARRAYLEAFSIZE / 2; ++i) {
		// split somehow:
		siblingNode -> keyArray[i] = oldLeafNode -> keyArray[(INTARRAYLEAFSIZE/2) + i]; // Right half of old node gets put into new node
		siblingNode -> ridArray[i] = oldLeafNode -> ridArray[(INTARRAYLEAFSIZE/2) + i];
		oldLeafNode -> keyArray[(INTARRAYLEAFSIZE/2) + i] = 0; // Right half of old node set to 0
		oldLeafNode -> ridArray[(INTARRAYLEAFSIZE/2) + i].page_number = 0;
	}

	// insert pair into the split leaves
	if (ridPair.key < siblingNode -> keyArray[0]) {
		insertLeaf(oldLeafNode, ridPair);
	}
	else {
		insertLeaf(siblingNode, ridPair);
	}

	// calculate new middle key pair:
	PageKeyPair<int> middleKeyPair;
	middleKeyPair.set(siblingId, siblingNode -> keyArray[0]);

	// return recursive to move up tree
	return middleKeyPair;
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
