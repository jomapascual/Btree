/**
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
	this->attrByteOffset = attrByteOffset;
	attributeType = attrType;
	scanExecuting = false;
	leafOccupancy = INTARRAYLEAFSIZE;
	nodeOccupancy = INTARRAYNONLEAFSIZE;

	std::ostringstream idxStr;
	idxStr << relationName << "." << attrByteOffset;
	outIndexName = idxStr.str();

	Page *headerPage, *rootPage;
	IndexMetaInfo* meta;

	try {
		file = new BlobFile(outIndexName, false); // Checks if the file exists

		headerPageNum = file -> getFirstPageNo();
		bufMgr -> readPage(file, headerPageNum, headerPage); // Reads the first page (meta)
		meta = (IndexMetaInfo*)headerPage;
		rootPageNum = meta -> rootPageNo; // Held in the meta page

		bufMgr -> unPinPage(file, headerPageNum, false); 
		//bufMgr -> unPinPage(file, rootPageNum, false);
	}
	catch (FileNotFoundException e) {
		file = new BlobFile(outIndexName, true);
		// Allocate header page
		bufMgr -> allocPage(file, headerPageNum, headerPage);
		memset(headerPage, 0, Page::SIZE);
		// Fills in meta info
		meta = (IndexMetaInfo*) headerPage;
		// Allocate root page
    	bufMgr -> allocPage(file, rootPageNum, (Page *&)rootPage); 
		memset(rootPage, 0, Page::SIZE);

		strcpy(meta -> relationName, relationName.c_str()); //strncpy
		meta -> attrByteOffset = attrByteOffset;
		meta -> attrType = attrType;
		meta -> rootPageNo = rootPageNum;
		initialRootPageNum = rootPageNum;
		
		LeafNodeInt *root = (LeafNodeInt *)rootPage; // Initializes the root
    	root -> rightSibPageNo = 0;

		bufMgr -> unPinPage(file, headerPageNum, true); 
		bufMgr -> unPinPage(file, rootPageNum, true); 

		// Scans Records
		FileScan scan = FileScan(relationName, bufMgr);
		//scanExecuting = true;
		
		try {
			RecordId rid;
			while(1) {
				scan.scanNext(rid);
				std::string fileRecord = scan.getRecord();
				const char *record = fileRecord.c_str();
				insertEntry(record + attrByteOffset, rid);
			}
		}
		catch (EndOfFileException e) {
			//bufMgr -> flushFile(file);
		}
	}
}


// -----------------------------------------------------------------------------
// BTreeIndex::~BTreeIndex -- destructor
// -----------------------------------------------------------------------------

BTreeIndex::~BTreeIndex()
{
	/// stop the ongoing scanning
	scanExecuting = false;
	bufMgr -> flushFile(BTreeIndex::file);
	delete file;
	file = nullptr;
}

/**
* Insert a key and a record ID to the appropriate position and array of the given leaf node.
*
* @param key   key to be inserted to keyArray
* @param rid   rid to be inserted to ridArray
* @param node  corresponding leaf node to be modified
* @param index index of a particular key and rid to be inserted
*/
void BTreeIndex::insertLeafNode(int key,
				    RecordId rid,
				    LeafNodeInt* node,
					int index) {
	// Shift to the right key after index and insert key at the appropriate position
	memmove(&node->keyArray[index + 1], &node->keyArray[index], sizeof(int) * (INTARRAYLEAFSIZE - index - 1));
	node->keyArray[index] = key;
	
	// Do the same for rid
	memmove(&node->ridArray[index + 1], &node->ridArray[index], sizeof(RecordId) * (INTARRAYLEAFSIZE - index - 1));
	node->ridArray[index] = rid;

	// Increment size
	node->size += 1;
}

/**
* Insert a key and a page ID to the appropriate position and array of the given non-leaf node.
*
* @param key   key to be inserted to keyArray
* @param childPageId   pageid to be inserted to pageNoArray
* @param node  corresponding non-leaf node to be modified
* @param index index of a particular key and pageid to be inserted
*/
void BTreeIndex::insertNonLeafNode(int key,
				    PageId childPageId,
					NonLeafNodeInt* node,
					int index) {
	// Shift to the right key after index and insert key at the appropriate position
	memmove(&node->keyArray[index + 1], &node->keyArray[index], sizeof(int) * (INTARRAYNONLEAFSIZE - index - 1));
	node->keyArray[index] = key;
	
	// Do the same for pid inside pageNoArray
	memmove(&node->pageNoArray[index + 2], &node->pageNoArray[index + 1], sizeof(PageId) * (INTARRAYNONLEAFSIZE - index - 1));
	node->pageNoArray[index + 1] = childPageId;

	// Increment size
	node->size += 1;
}

/**
* Helper function to do the recurion of inserting a key and a rid to the tree.
*
* @param key   key to be inserted
* @param rid   rid to be inserted
* @param currPageId  the current page which holds the node
* @param middleValueFromChild a callback parameter to pass back the middle value to its parent
* @param newlyCreatedPageId a callback parameter to pass back the splitted node's pageID to its parent
* @param isLeafBool whether this node is a leaf or not
*/
void BTreeIndex::insertEntryHelper(const int key,
					const RecordId rid,
					PageId currPageId,
					int* middleValueFromChild,
					PageId* newlyCreatedPageId,
					bool isLeafBool) 
{
	// Read page, which is a node
	Page *currNode;
  	bufMgr->readPage(file, currPageId, currNode);

  	// Base Case: current node is a leaf node
  	if (isLeafBool) {
  		LeafNodeInt* currLeafNode = (LeafNodeInt *)currNode;

  		// Index of the new key to be inserted
  		int index = currLeafNode->size;
  		for(int i = 0; i < currLeafNode->size; i++) {
  			if (currLeafNode->keyArray[i] > key) {
  				index = i;
  				break;
  			}
  		}

  		// Check if leaf is not full, if true directly insert to the leaf
  		if (currLeafNode->ridArray[INTARRAYLEAFSIZE - 1].page_number == Page::INVALID_NUMBER) {
  		//if (currLeafNode->size < INTARRAYLEAFSIZE) {	
			insertLeafNode(key, rid, currLeafNode, index);
  			bufMgr->unPinPage(file, currPageId, true);
  			*middleValueFromChild = 0;
  			*newlyCreatedPageId = 0;
  		} 
  		// Otherwise, split node into 2 and pass back the middle value
  		else {
  			splitLeaf(currLeafNode, currPageId, key, rid, index, middleValueFromChild, newlyCreatedPageId);
		}
  	}
  	// Recursive Case: current node is not a leaf node
  	else {
  		NonLeafNodeInt *currNonLeafNode = (NonLeafNodeInt *)currNode;

  		// Find the correct child node
  		int childIndex = currNonLeafNode->size;
  		for(int i = 0; i < currNonLeafNode->size; i++) {
  			if (currNonLeafNode->keyArray[i] > key) {
  				childIndex = i;
  				break;
  			}
  		}
  		PageId currChildId = currNonLeafNode->pageNoArray[childIndex];

  		// Recursive call to the child
		int newChildMiddleKey;
		PageId newChildId;
		insertEntryHelper(key, rid, currChildId, &newChildMiddleKey, &newChildId, currNonLeafNode->level == 1);

		// If there is no split in child node
		if ((int) newChildId == 0) {
		  bufMgr->unPinPage(file, currPageId, false);
		  *middleValueFromChild = 0;
		  *newlyCreatedPageId = 0;
		}
		// If there is a split in child node
		else {
	  		// Index of the new middle key from children to be inserted
	  		int index = currNonLeafNode->size;
	  		for(int i = 0; i < currNonLeafNode->size; i++) {
	  			if (currNonLeafNode->keyArray[i] > newChildMiddleKey) {
	  				index = i;
	  				break;
	  			}
	  		}

	  		// Check if node is not full, if true directly insert middle key to the node
  			if (currNonLeafNode->pageNoArray[INTARRAYNONLEAFSIZE] == Page::INVALID_NUMBER) {
  			//if (currNonLeafNode->size < INTARRAYNONLEAFSIZE) {
  				insertNonLeafNode(newChildMiddleKey, newChildId, currNonLeafNode, index);
  				bufMgr->unPinPage(file, currPageId, true);
	  			*middleValueFromChild = 0;
	  			*newlyCreatedPageId = 0;
  			}
  			// Otherwise, split node into 2 and pass back the middle value
  			else {
	  			// Allocate a new page for the right split
	  			PageId newPageId;
	  			Page* newNode;
	  			bufMgr->allocPage(file, newPageId, newNode);
	  			memset(newNode, 0, Page::SIZE);
	  			NonLeafNodeInt* newNonLeafNode = (NonLeafNodeInt *)newNode;
	  			newNonLeafNode->level = currNonLeafNode->level;

	  			// Split nodes depending on the cases
	  			int mid = INTARRAYNONLEAFSIZE / 2;
	  			// Case 1: the new child will be pushed
	  			if (index == mid) {
	  				for(int i = mid; i < INTARRAYNONLEAFSIZE; i++) {
						newNonLeafNode->keyArray[i-mid] = currNonLeafNode->keyArray[i];
						newNonLeafNode->pageNoArray[i-mid+1] = currNonLeafNode->pageNoArray[i+1];
						currNonLeafNode->keyArray[i] = 0;
						currNonLeafNode->pageNoArray[i+1] = Page::INVALID_NUMBER;
					}
					newNonLeafNode->pageNoArray[0] = newChildId;

					// Set size appropriately
					currNonLeafNode->size = mid;
					newNonLeafNode->size = INTARRAYNONLEAFSIZE - mid;

		  			// Unpin the nodes
		  			bufMgr->unPinPage(file, currPageId, true);
		  			bufMgr->unPinPage(file, newPageId, true);

					// Return values using pointers
		  			*middleValueFromChild = newChildMiddleKey;
		  			*newlyCreatedPageId = newPageId;
	  			}
	  			// Case 2: if new child will not be pushed
	  			else {
	  				if (INTARRAYNONLEAFSIZE % 2 == 0 && index < mid) {
	  					mid -= 1;
	  				}
	  				for(int i = mid + 1; i < INTARRAYNONLEAFSIZE; i++) {
						newNonLeafNode->keyArray[i-mid-1] = currNonLeafNode->keyArray[i];
						newNonLeafNode->pageNoArray[i-mid-1] = currNonLeafNode->pageNoArray[i];
						currNonLeafNode->keyArray[i] = 0;
						currNonLeafNode->pageNoArray[i] = Page::INVALID_NUMBER;
					}

					// Return values using pointers
					*middleValueFromChild = currNonLeafNode->keyArray[mid];
					*newlyCreatedPageId = newPageId;

					// Clear already pushed value
					currNonLeafNode->keyArray[mid] = 0;

					// Set size appropriately
					currNonLeafNode->size = mid;
					newNonLeafNode->size = INTARRAYNONLEAFSIZE - mid - 1;

					// Insert new value to left or right node appropriately
					if (index < INTARRAYNONLEAFSIZE / 2) {
						insertNonLeafNode(newChildMiddleKey, newChildId, currNonLeafNode, index);
					} else {
						insertNonLeafNode(newChildMiddleKey, newChildId, newNonLeafNode, index - mid);
					}

		  			// Unpin the nodes
		  			bufMgr->unPinPage(file, currPageId, true);
		  			bufMgr->unPinPage(file, newPageId, true);
	  			}
  			}
		}
  	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
{
	// Call the helper function to do the recursion while retrieving back new middle value and pageId if there is a split
	int middleValueFromChild;
	PageId newlyCreatedPageId;
	insertEntryHelper(*(int*)key, rid, rootPageNum, &middleValueFromChild, &newlyCreatedPageId, rootIsLeaf);

	// If there is a split to the root, create a new root
	if ((int) newlyCreatedPageId != 0) {
	  	// Allocate a new page for this new root
	  	PageId newPageId;
		Page* newPage;
		bufMgr->allocPage(file, newPageId, newPage);
		memset(newPage, 0, Page::SIZE);
		NonLeafNodeInt* newRoot = (NonLeafNodeInt *)newPage;
		
		// Update the new page appropriately
		newRoot->keyArray[0] = middleValueFromChild;
		newRoot->pageNoArray[0] = rootPageNum;
		newRoot->pageNoArray[1] = newlyCreatedPageId;
		newRoot->size = 1;
		if(rootIsLeaf) {
			newRoot->level = 1;
		} else {
			newRoot->level = 0;
		}
		rootIsLeaf = false;

		// Update global variable and IndexMetaInfo page appropriately
		Page *meta;
		bufMgr->readPage(file, headerPageNum, meta);
		IndexMetaInfo *metadata = (IndexMetaInfo *)meta;
		metadata->rootPageNo = newPageId;
		rootPageNum = newPageId;

		// Unpin the root and the IndexMetaInfo page
		bufMgr->unPinPage(file, newPageId, true);
		bufMgr->unPinPage(file, headerPageNum, true);
	}
}

const void BTreeIndex::splitLeaf(LeafNodeInt *currLeafNode, PageId leafPageNum, const int key, const RecordId rid, const int index, int* middleValueFromChild, PageId* newlyCreatedPageId){
	PageId newPageId;
	Page* newNode;
	bufMgr->allocPage(file, newPageId, newNode);
	memset(newNode, 0, Page::SIZE);
	LeafNodeInt* newLeafNode = (LeafNodeInt *)newNode;

	// Split and add new value appropriately depending on the position of the index
	int mid = INTARRAYLEAFSIZE / 2;
	if (INTARRAYLEAFSIZE % 2 == 1 && index > mid) {
		mid = mid + 1;
	}

	for(int i = mid; i < INTARRAYLEAFSIZE; i++) {
		newLeafNode->keyArray[i-mid] = currLeafNode->keyArray[i];
		newLeafNode->ridArray[i-mid] = currLeafNode->ridArray[i];
		currLeafNode->keyArray[i] = 0;
		currLeafNode->ridArray[i].page_number = Page::INVALID_NUMBER;
	}

	// Set size appropriately
	newLeafNode->size = INTARRAYLEAFSIZE - mid;
	currLeafNode->size = mid;

	// Insert to right node
	if(index > INTARRAYLEAFSIZE / 2) {
		insertLeafNode(key, rid, newLeafNode, index - mid);
	}
	// Insert to left node
	else {
		insertLeafNode(key, rid, currLeafNode, index);
	}

	// Set rightSibPageNo appropriately
	newLeafNode->rightSibPageNo = currLeafNode->rightSibPageNo;
	currLeafNode->rightSibPageNo = newPageId;

	// Unpin the nodes
	bufMgr->unPinPage(file, leafPageNum, true);
	bufMgr->unPinPage(file, newPageId, true);

	*middleValueFromChild = newLeafNode->keyArray[0];
	*newlyCreatedPageId = newPageId;	
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
const int BTreeIndex::startScanHelper(PageId pageNum)
{
	NonLeafNodeInt* currNode;
	LeafNodeInt* child;

	currentPageNum = pageNum;
	bufMgr -> readPage(file, currentPageNum, currentPageData);
	currNode = (NonLeafNodeInt*)currentPageData;
	for (int i = 0; i < INTARRAYNONLEAFSIZE + 1; ++i) {
		if (i == INTARRAYNONLEAFSIZE || currNode -> pageNoArray[i + 1] == Page::INVALID_NUMBER || currNode -> keyArray[i] > lowValInt) {
            // node is directly above leaf node 
            if (currNode -> level == 1) { 
                // read child page (leaf), then update
                currentPageNum = currNode -> pageNoArray[i];
                bufMgr -> readPage(file, currentPageNum, currentPageData);
                child = (LeafNodeInt*)currentPageData;
                // scan page, return index
                for(int j = 0; j < INTARRAYLEAFSIZE; ++j) {
                    if (child -> ridArray[j].page_number == Page::INVALID_NUMBER) {
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
	int nextIndex;
	if(leafNode->ridArray[nextEntry].page_number == Page::INVALID_NUMBER) { // Checks for valid nextEntry
		throw IndexScanCompletedException();
	}
	else if(nextEntry == INT_MAX) {
		throw IndexScanCompletedException(); // Scan has been completed
	}
	else {
		nextIndex = nextEntry;
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
		if(leafNode -> ridArray[nextEntry + 1].page_number != Page::INVALID_NUMBER) { // Checks if there are no more pages left
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