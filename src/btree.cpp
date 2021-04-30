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
/**
 * BTreeIndex Constructor. 
 * Check to see if the corresponding index file exists. If so, open the file.
 * If not, create it and insert entries for every tuple in the base relation using FileScan class.
 *
 * @param relationName      Name of file.
 * @param outIndexName      Return the name of index file.
 * @param bufMgrIn		Buffer Manager Instance
 * @param attrByteOffset	Offset of attribute, over which index is to be built, in the record
 * @param attrType		Datatype of attribute over which index is built
 * @throws  BadIndexInfoException     If the index file already exists for the corresponding attribute, but values in metapage(relationName, attribute byte offset, attribute type etc.) do not match with values received through constructor parameters.
 */
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
	idxStr << relationName << "." << attrByteOffset;
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
		initialRootPageNum = rootPageNum;
		
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
/**
 * BTreeIndex Destructor. 
 * End any initialized scan, flush index file, after unpinning any pinned pages, from the buffer manager
 * and delete file instance thereby closing the index file.
 * Destructor should not throw any exceptions. All exceptions should be caught in here itself. 
 * @param none
 * @return void
 */
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
	// RIDKeyPair<const void*> newEntry;
	// newEntry.set(rid, key);
	// Page *currPage;
	// NonLeafNodeInt* nonLeafNode = (LeafNodeInt *)currPage;
	// bool split = false;

	// // leaf
	// if (nonLeafNode -> level = -1) {
	// 	bufMgr -> unPinPage(file, pid, false);
	// 	split = insertLeaf(currPage, newEntry); 	// not sure
	// }
	// // nonLeaf
	// else {
	// 	// find index
	// 	int i = 0;
	// 	while(node->pageNoArray[i+1] != 0 && key > node->keyArray[i] && i < INTARRAYNONLEAFSIZE) {
	// 		i++;
	// 	}
	// 	index = i;
	// 	PageId childPid = nonLeaf->pageNoArray[index];
	// 	bufMgr->unPinPage(currPage, pageNum, true);
	// 	//recursively find the leaf node, insert a pushed up entry to current node if son is splitted
		
	// }
	RIDKeyPair<int> newEntry;
	newEntry.set(rid, *((int *)key);
	Page* rootPage;
	bufMgr -> readPage(file, rootPageNum, rootPage);
  	PageKeyPair<int> *newChild;
	
	bool isLeaf;
	if(initialRootPageNum == rootPageNum) {
		isLeaf = true;
	}
	else {
		isLeaf= false;
	}
	insertHelper(rootPage, rootPageNum, isLeaf, newEntry, newChild);

}

void insertHelper(Page *page, PageId pageNum, bool isLeaf, const RIDKeyPair<int> newEntry, PageKeyPair<int> *&child) {
	if(initialRootPageNum == rootPageNum) { // If node is leaf
		LeafNodeInt *leaf = (LeafNodeInt *)page;
    
		if (leaf -> ridArray[INTARRAYLEAFSIZE - 1].page_number == 0)	{ // Page is not full
			insertLeaf(leaf, newEntry);
			bufMgr -> unPinPage(file, pageNum, true);
			child = nullptr;
		}
		else {
			splitLeafNode(leaf, pageNum, newEntry.key); // might need PageKeyPair
		}
	}
	else {	// If node is not leaf
		NonLeafNodeInt *currNode = (NonLeafNodeInt *)page;
		Page *nextPage;
		PageId nextNodeIndex;

		int i = INTARRAYNONLEAFSIZE; // Will be the index of the next nonLeafNode
		while ((currNode -> keyArray[i-1] >= child -> key) && i > 0) { 
			i--;
		}
		while ((currNode -> pageNoArray[i] == 0) && i >= 0) {
			i--;
		}
		nextNodeIndex = currNode -> pageNoArray[i];
		bufMgr -> readPage(file, nextNodeIndex, nextPage);

		isLeaf;
		if(currNode -> level == 1){
			isLeaf = true;
		}
		else {
			isLeaf = false;
		}
		insertHelper(nextPage, pageNum, isLeaf, newEntry, child);

		if(child != nullptr) { // Child has a split
			if (currNode -> pageNoArray[INTARRAYNONLEAFSIZE] != 0) { // Current page is full
				splitNonLeafNode(currNode, pageNum, child);
			}
			else {
				insertNonLeafNode(currNode, child); // Inserts PageKeypair to the current page
				bufMgr -> unPinPage(file, pageNum, true); // Unpins reaching the end of insert
				child = nullptr;
			}
		}
	}	
}

/**
 * A helper method for inserting a leaf node in the B Tree
 * @param leafNode pointer to the leaf node to insert
 * @param ridKey the RIDKeyPair of the node to insert
 */
void insertLeaf(LeafNodeInt *leafNode, RIDKeyPair<int> ridKey) // Should use parameters like this
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

/**
 * A helper method for inserting a nonleaf node in the B Tree
 * @param node pointer to the nonleaf node to insert
 * @param KeyPage the PageKeyPair of the node to insert
 */
void insertNonLeaf(NonLeafNodeInt *node, PageKeyPair<int> keyPage)
{
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
	
/**
 * A helper method for splitting a leaf node in the B Tree
 * @param oldLeafNode pointer to the leaf node that needs to be split
 * @param currPageId the PageID of the current leaf node
 * @param ridPair the RIDKeyPair of the given leaf node
 * @return PageKeyPair<int> 
 */
void splitLeafNode(LeafNodeInt *oldLeafNode, PageId currPageId, RIDKeyPair<int> ridPair) //FIXME don't need to return anything
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

	siblingNode -> rightSibPageNo = oldLeafNode -> rightSibPageNo; // Sets the sibling pointer
	oldLeafNode -> rightSibPageNo = siblingId;

	// calculate new middle key pair:
	PageKeyPair<int> middleKeyPair;
	middleKeyPair.set(siblingId, siblingNode -> keyArray[0]);

	// the smallest key from second page as the new child entry
	middleKeyPair = new PageKeyPair<int>();
	PageKeyPair<int> newKeyPair;
	newKeyPair.set(siblingId, newLeafNode->keyArray[0]);
	middleKeyPair = &newKeyPair;
	bufMgr->unPinPage(file, currPageId, true);
	bufMgr->unPinPage(file, siblingId, true);

	// if curr page is root
	if (currPageId == rootPageNum) {
		updateRoot(currPageId, middleKeyPair);
	}
  }

/**
 * A helper method for splitting a nonleaf node in the B Tree
 * @param oldNonLeafNode pointer to the nonleaf node that needs to be split
 * @param currPageId the PageID of the current nonleaf node
 * @param newPair the PageKeyPair of the given nonleaf node
 */
void splitNonLeafNode(NonLeafNodeInt *oldNonLeafNode, PageId currPageId, PageKeyPair<int>*&newPair)
 {
	 // allocate new non-leaf sibling
	 Page* sibling;
	 PageId siblingId;
	 bufMgr -> allocPage(file, siblingId, sibling);
	 // convert to proper structure
	 NonLeafNodeInt* siblingNode = (NonLeafNodeInt*) sibling;

	 // middle range
	 int middle = INTARRAYNONLEAFSIZE / 2;
	 int pushUp = middle;
	 PageKeyPair<int> nodeToPush;

	// check if full, to push up
	 if (INTARRAYNONLEAFSIZE % 2 == 0) {
		 if (newPair -> key < oldNonLeafNode -> keyArray[middle]) {
			 pushUp = middle - 1;
		 } else {
			 pushUp = middle;
		 }
	 }
	 nodeToPush.set(siblingId, oldNonLeafNode->keyArray[pushUp]);

	middle = pushUp + 1;
	for (int i = middle; i < INTARRAYNONLEAFSIZE; i++) {
		siblingNode -> keyArray[i - middle] = oldNonLeafNode -> keyArray[i];
		siblingNode -> pageNoArray[i - middle] = oldNonLeafNode -> pageNoArray[i + 1];
		oldNonLeafNode -> pageNoArray[i + 1] = (PageId) 0;
		oldNonLeafNode -> keyArray[i + 1] = 0;
	}

	siblingNode -> level = oldNonLeafNode -> level;
	oldNonLeafNode -> keyArray[pushUp] = 0;
	oldNonLeafNode -> pageNoArray[pushUp] = (PageId) 0;

	if (newPair -> key < siblingNode -> keyArray[0]) {
		insertNonLeaf(oldNonLeafNode, *newPair);
	} 
	else { 
		insertNonLeaf(siblingNode, *newPair);
	}

	newPair = &nodeToPush;
	bufMgr -> unPinPage(file, currPageId, true);
	bufMgr -> unPinPage(file, siblingId, true);

	if (currPageId == rootPageNum) {
		updateRootNode(*newPair, currPageId);
	}
 }

/**
 * A helper method for updating the root in the B Tree
 * @param keyPage the PageKeyPair of the root node
 * @param currPage the PageID of the given node
 */
void updateRootNode(PageKeyPair<int> keyPage, PageId currPage) {
	Page* root;
	PageId rootId; // new root page number
	bufMgr -> allocPage(file, rootId, root);
	// convert to proper structure
	NonLeafNodeInt *rootPage = (NonLeafNodeInt*) root;

	// Need update level ?
	rootPage -> keyArray[0] = keyPage.key;
	rootPage -> pageNoArray[0] = currPage;
	rootPage -> pageNoArray[1] = keyPage.pageNo;

	Page *headerInfo;
	bufMgr -> readPage(file, headerPageNum, headerInfo); // Reads the first page (meta)
	IndexMetaInfo *headerPage = (IndexMetaInfo *) headerInfo;
	rootPageNum = rootId;
	headerPage -> rootPageNo = rootId;

	// Update pins needed?
	bufMgr -> unPinPage(file, rootId, true);
}


// -----------------------------------------------------------------------------
// BTreeIndex::startScan
// -----------------------------------------------------------------------------
/**
 * Begin a filtered scan of the index.  For instance, if the method is called 
 * using ("a",GT,"d",LTE) then we should seek all entries with a value 
 * greater than "a" and less than or equal to "d".
 * If another scan is already executing, that needs to be ended here.
 * Set up all the variables for scan. Start from root to find out the leaf page that contains the first RecordID
 * that satisfies the scan parameters. Keep that page pinned in the buffer pool.
 * @param lowVal	Low value of range, pointer to integer / double / char string
 * @param lowOp		Low operator (GT/GTE)
 * @param highVal	High value of range, pointer to integer / double / char string
 * @param highOp	High operator (LT/LTE)
 * @throws  BadOpcodesException If lowOp and highOp do not contain one of their their expected values 
 * @throws  BadScanrangeException If lowVal > highval
 * @throws  NoSuchKeyFoundException If there is no key in the B+ tree that satisfies the scan criteria.
 */
void BTreeIndex::startScan(const void* lowValParm,
				   const Operator lowOpParm,
				   const void* highValParm,
				   const Operator highOpParm)
{
    // Add your code below. Please do not remove this line.

	// end other scans
	if (scanExecuting) {
		endScan();
	}

	// initial conditions and exceptions
	if (lowOpParm != GT && lowOpParm != GTE) {
		throw BadOpcodesException();
	}
	if (highOpParm != LT && highOpParm != LTE) {
		throw BadOpcodesException();
	}
	
	lowValInt = *((int *)lowValParm);
	highValInt = *((int *)highValParm);
	if (lowValInt > highValInt) {
		throw BadScanrangeException();
	}

	lowOp = lowOpParm;
	highOp = highOpParm;

	// start new scan
	scanExecuting = true;
	nextEntry = startScanHelper(rootPageNum);

}

/**
 * Helper method for the startScan method to update nextEntry number
 * @param pageNum 	 
 * @return 	int 
 */
const int startScanHelper(PageId pageNum)
{
	NonLeafNodeInt* currNode;
	LeafNodeInt* child;

	currentPageNum = pageNum;
	bufMgr -> readPage(file, currentPageNum, currentPageData);
	currNode = (NonLeafNodeInt*)currentPageData;
	for (int i = 0; i < INTARRAYNONLEAFSIZE + 1; ++i) {
		if (i == INTARRAYNONLEAFSIZE || currNode -> pageNoArray[i + 1] == PAGE::INVALID_NUMBER || currNode -> keyArray[i] > lowValInt) {
            // node is directly above leaf node 
            if (currNode -> level == 1) { 
                // read child page (leaf), then update
                currentPageNum = currNode -> pageNoArray[i];
                bufMgr -> readPage(file, currentPageNum, currentPageData);
                child = (LeafNodeInt*)currentPageData;
                // scan page, return index
                for(int j = 0; j < INTARRAYLEAFSIZE; ++j) {
                    if (childNode -> ridArray[j] == Page::INVALID_SLOT) {
                        break;
                    }
                    if ((lowOp == GT && child -> keyArray[j] > lowValInt)
                    || (lowOp == GTE && child -> keyArray[j] >= lowValInt)) {
                        bufMgr -> unPinPage(file, pageNum, false);
                        return j;
                    }
                }   
                // if not found by the end
                bufMgr -> unPinPage(file, pageNum, false); 
                endScan();  
                throw NoSuchKeyFoundException();      
            } 
			else { 
				// recursive call for nonLeaf nodes
                try {
                    int next = startScanHelper(currNode -> pageNoArray[i]);
                    bufMgr -> unPinPage(file, pageNum, false);
                    return next;
                } catch (NoSuchKeyFoundException e) {
                    bufMgr -> unPinPage(file, pageNum, false);
                    throw NoSuchKeyFoundException(); 
                }
            }
        }
    }

	// impossible
    assert(false);
    return -1;

}


// -----------------------------------------------------------------------------
// BTreeIndex::scanNext
// -----------------------------------------------------------------------------
/**
 * Fetch the record id of the next index entry that matches the scan.
 * Return the next record from current page being scanned. If current page has been scanned to its entirety, 
 * move on to the right sibling of current page, if any exists, to start scanning that page. Make sure to unpin any pages that are no longer required.
 * @param outRid	RecordId of next record found that satisfies the scan criteria returned in this
 * @throws ScanNotInitializedException If no scan has been initialized.
 * @throws IndexScanCompletedException If no more records, satisfying the scan criteria, are left to be scanned.
 **/
void BTreeIndex::scanNext(RecordId& outRid) 
{
    // Add your code below. Please do not remove this line.
	if(scanExecuting == false) { // Checks if scan has been initialized
		throw ScanNotInitializedException();
	}

	LeafNodeInt* leafNode = (LeafNodeInt*) currentPageData;
	if(leafNode->ridArray[nextEntry] == INVALID_RECORD) { // Checks for valid nextEntry
		throw IndexScanCompletedException();
	}
	else if(nextEntry == INT_MAX) {
		throw IndexScanCompletedException(); // Scan has been completed
	}
	else {
		int nextIndex = nextEntry;
	}

	if((leafNode -> keyArray[nextIndex] <= highValInt) && highOp == LTE) {
		outRid = leafNode -> ridArray[nextIndex];
	}
	else if((leafNode -> keyArray[nextIndex] < highValInt) && highOp == LT) {
		outRid = leafNode -> ridArray[nextIndex];
	}
	else {
		throw IndexScanCompletedException(); // Scan has been completed
	}

	if(nextEntry < INTARRAYLEAFSIZE - 1) {
		if(leafNode -> ridArray[nextEntry + 1] != INVALID_RECORD) { // Checks if there are no more pages left
			nextEntry++; // increments nextEntry index
		}
		else {
			if(leafNode -> rightSibPageNo == Page::INVALID_NUMBER) {
				nextEntry = INT_MAX;
			}
			else { // Page has been fully scanned
				PageId nextPageId = currentPageNum;
				currentPageNum = leafNode -> rightSibPageNo; // Advances to the right sibling node
				nextEntry = 0;
				bufMgr -> readPage(file, currentPageNum, currentPageData); 
				bufMgr -> unPinPage(file, nextPageId, false); // Unpins fully scanned pages
			}
		}
	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::endScan
// -----------------------------------------------------------------------------
//
/**
 * Terminate the current scan. Unpin any pinned pages. Reset scan specific variables.
 * @throws ScanNotInitializedException If no scan has been initialized.
**/
void BTreeIndex::endScan() 
{
    // Add your code below. Please do not remove this line.
	//check if scan has been initialized
	if(scanExecuting == false){
		throw ScanNotInitializedException();
	}
	//terminate scan
	scanExecuting = false;
	//unpin pages
	bufMgr -> unPinPage(file, currentPageNum, false);
	//reset scan specific variables
	nextEntry = -1;
	currentPageNum = -1;
	currentPageData = nullptr;
	lowValInt = -1;
	highValInt = -1;
	lowOp = (Operator)-1;
	highOp = (Operator)-1;
}

}
